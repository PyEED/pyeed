from __future__ import annotations

import asyncio
from collections.abc import AsyncIterator
from time import monotonic
from typing import Final, Literal

import numpy as np
from loguru import logger
from neo4j import AsyncDriver
from pymilvus import AsyncMilvusClient
from rich.progress import Progress, TaskID

from ..db.milvus import (
    get_async_milvus_client,
    initialize_collection_from_dict,
    insert,
)
from ..db.neo4j import get_async_driver
from .new_embed import DType, ESM2Processor

__all__ = [
    "run_embedding_job",
    "stream_pending_proteins",
    "update_protein_status",
]

type ProteinRow = tuple[str, str]
type VectorResult = tuple[str, np.ndarray]
type StatusLiteral = Literal["pending", "in_progress", "complete", "failed"]

_PENDING_QUERY: Final[str] = """
MATCH (p:Protein {embedding_status: 'pending'})
RETURN p.id AS id, p.sequence AS sequence
"""


async def update_protein_status(
    driver: AsyncDriver,
    ids: list[str],
    status: StatusLiteral,
) -> None:
    """Update embedding_status for a list of proteins."""
    if not ids:
        return

    query = """
    UNWIND $ids AS id
    MATCH (p:Protein {id: id})
    SET p.embedding_status = $status
    """
    async with driver.session() as session:
        await session.run(query, ids=ids, status=status)

    logger.debug(
        "Updated protein statuses",
        extra={"count": len(ids), "status": status},
    )


async def stream_pending_proteins(
    driver: AsyncDriver,
    *,
    batch_size: int = 5_000,
    mark_in_progress: bool = True,
) -> AsyncIterator[list[ProteinRow]]:
    """Stream pending protein ids and sequences in bounded batches.

    Keeps the Neo4j session open for streaming to avoid loading all rows into
    memory while still yielding manageable batches for the embedder.
    Marks yielded proteins as ``in_progress`` when mark_in_progress is True.
    """
    async with driver.session() as session:
        result = await session.run(_PENDING_QUERY)

        buffer: list[ProteinRow] = []
        async for record in result:
            protein_id = record["id"]
            sequence = record["sequence"]
            buffer.append((protein_id, sequence))

            if len(buffer) >= batch_size:
                logger.debug(
                    "Yielding pending batch",
                    extra={"batch_size": len(buffer)},
                )
                if mark_in_progress:
                    await update_protein_status(driver, [pid for pid, _ in buffer], "in_progress")
                yield buffer
                buffer = []

        if buffer:
            logger.debug(
                "Yielding final pending batch",
                extra={"batch_size": len(buffer)},
            )
            if mark_in_progress:
                await update_protein_status(driver, [pid for pid, _ in buffer], "in_progress")
            yield buffer


async def run_embedding_job(
    *,
    neo4j_driver: AsyncDriver,
    milvus_client: AsyncMilvusClient,
    collection_name: str,
    devices: list[int] | None = None,
    dtype: DType = "float16",
    max_padded_tokens: int = 30_000,
    max_length: int = 1_024,
    source_chunk_size: int = 2_048,
    read_batch_size: int = 5_000,
    writer_queue_size: int = 2_000,
    writer_batch_size: int = 1_000,
    status_batch_size: int = 1_000,
    huggingface_token: str | None = None,
    progress: Progress | None = None,
    progress_task_id: TaskID | None = None,
) -> None:
    """Run the end-to-end embedding job with backpressure and streaming IO.

    Streams pending proteins from Neo4j, marks them in progress, embeds with
    ESM2, writes embeddings to Milvus, and updates statuses to complete/failed.
    """
    out_q: asyncio.Queue[VectorResult | None] = asyncio.Queue(writer_queue_size)
    collection_initialized = collection_name in await milvus_client.list_collections()

    async def flush(batch: list[VectorResult]) -> None:
        nonlocal collection_initialized
        if not batch:
            return
        batch_count = len(batch)
        batch_ids = [pid for pid, _ in batch]
        start = monotonic()
        records = [{"id": pid, "embedding": emb} for pid, emb in batch]

        try:
            if not collection_initialized:
                await initialize_collection_from_dict(
                    collection_name=collection_name,
                    primary_field_name="id",
                    record=records[0],
                    client=milvus_client,
                )
                collection_initialized = True
                logger.info(
                    "Initialized Milvus collection",
                    extra={"collection": collection_name},
                )

            await insert(
                collection_name=collection_name,
                records=records,
                client=milvus_client,
            )
            await update_protein_status(neo4j_driver, batch_ids, "complete")

            logger.debug(
                "Writer flushed batch",
                extra={
                    "items": batch_count,
                    "elapsed_s": round(monotonic() - start, 3),
                },
            )
        except Exception:
            logger.exception(
                "Failed to flush batch to Milvus",
                extra={"items": batch_count},
            )
            try:
                await update_protein_status(neo4j_driver, batch_ids, "failed")
            except Exception:
                logger.exception(
                    "Failed to mark proteins as failed after Milvus error",
                    extra={"count": len(batch_ids)},
                )
        finally:
            batch.clear()

    async def writer() -> None:
        buffer: list[VectorResult] = []
        backlog_warned = False
        while True:
            item = await out_q.get()
            try:
                if item is None:
                    break
                buffer.append(item)

                qsize = out_q.qsize()
                maxsize = out_q.maxsize
                if maxsize and qsize > maxsize * 0.8 and not backlog_warned:
                    backlog_warned = True
                    logger.debug(
                        "Writer backlog high",
                        extra={"qsize": qsize, "maxsize": maxsize},
                    )
                elif backlog_warned and maxsize and qsize < maxsize * 0.5:
                    backlog_warned = False

                if len(buffer) >= writer_batch_size:
                    await flush(buffer)
                    # batch cleared inside flush
                if len(buffer) >= status_batch_size:
                    await flush(buffer)
            finally:
                out_q.task_done()

        if buffer:
            await flush(buffer)

    chunks = stream_pending_proteins(neo4j_driver, batch_size=read_batch_size)
    writer_task = asyncio.create_task(writer())

    async with ESM2Processor(
        devices=devices,
        max_padded_tokens=max_padded_tokens,
        max_length=max_length,
        source_chunk_size=source_chunk_size,
        dtype=dtype,
        huggingface_token=huggingface_token,
        progress=progress,
        progress_task_id=progress_task_id,
    ) as processor:
        async for protein_id, embedding in processor.stream_work(chunks):
            put_start = monotonic()
            await out_q.put((protein_id, embedding))
            waited = monotonic() - put_start
            if waited > 0.1:
                logger.debug(
                    "Backpressure on embedding queue",
                    extra={
                        "wait_s": round(waited, 3),
                        "qsize": out_q.qsize(),
                        "queue_max": out_q.maxsize,
                    },
                )

    await out_q.put(None)
    await out_q.join()
    await writer_task


async def main() -> None:
    """Example entry point for running the embedding job."""
    driver = get_async_driver()
    milvus_client = get_async_milvus_client()
    try:
        # Progress total is unknown until we count; use indeterminate progress.
        with Progress() as progress:
            task_id = progress.add_task("embedding", total=None)
            await run_embedding_job(
                neo4j_driver=driver,
                milvus_client=milvus_client,
                collection_name="proteins",
                progress=progress,
                progress_task_id=task_id,
                devices=[0, 2],
            )
    finally:
        await milvus_client.close()
        await driver.close()


if __name__ == "__main__":
    asyncio.run(main())
