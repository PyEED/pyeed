from __future__ import annotations

import argparse
import asyncio
from collections.abc import AsyncIterator
from time import monotonic
from typing import Final, Literal

import numpy as np
from loguru import logger
from neo4j import AsyncDriver
from pymilvus import AsyncMilvusClient
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskID,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

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

_STATUS_COUNTS_QUERY: Final[str] = """
MATCH (p:Protein)
RETURN p.embedding_status AS status, count(*) AS count
"""


async def count_pending(driver: AsyncDriver) -> tuple[int, dict[str, int]]:
    """Return pending count and counts by status for logging/progress."""
    async with driver.session() as session:
        result = await session.run(_STATUS_COUNTS_QUERY)
        rows = await result.data()
    counts = {row["status"]: row["count"] for row in rows}
    return int(counts.get("pending", 0)), counts


def parse_devices(devices_str: str) -> list[int]:
    """Parse comma-separated CUDA device IDs, or 'cpu' for CPU-only."""
    if not devices_str:
        return []
    lowered = devices_str.strip().lower()
    if lowered == "cpu":
        return [-1]  # sentinel to force CPU in ESM2Processor
    try:
        return [int(d.strip()) for d in devices_str.split(",") if d.strip()]
    except ValueError as exc:
        msg = "Invalid --devices value. Use comma-separated integers like '0,1' or 'cpu'."
        raise ValueError(msg) from exc


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
    pending_total, status_counts = await count_pending(neo4j_driver)
    logger.info(
        "Embedding job starting",
        extra={"pending": pending_total, "status_counts": status_counts},
    )
    if pending_total == 0:
        logger.info("No pending proteins to embed; exiting")
        return
    if progress is not None and progress_task_id is not None:
        progress.update(progress_task_id, total=pending_total)

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


async def main(argv: list[str] | None = None) -> None:
    """CLI entry point for running the embedding job."""
    parser = argparse.ArgumentParser(description="Embed pending proteins from Neo4j into Milvus")
    parser.add_argument(
        "--collection",
        type=str,
        default="proteins",
        help="Milvus collection name (default: proteins)",
    )
    parser.add_argument(
        "--devices",
        type=str,
        default="0",
        help="Comma-separated CUDA device IDs (e.g. '0,1') or 'cpu' (default: 0)",
    )
    parser.add_argument(
        "--dtype",
        type=str,
        choices=["float16", "float32"],
        default="float16",
        help="Embedding dtype (default: float16)",
    )
    parser.add_argument(
        "--max-padded-tokens",
        type=int,
        default=30_000,
        help="Max padded tokens per batch for ESM2 (default: 30000)",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=1_024,
        help="Max tokens per sequence (default: 1024)",
    )
    parser.add_argument(
        "--source-chunk-size",
        type=int,
        default=2_048,
        help="Chunk size for source stream into the embedder (default: 2048)",
    )
    parser.add_argument(
        "--read-batch-size",
        type=int,
        default=5_000,
        help="Rows to read from Neo4j per batch (default: 5000)",
    )
    parser.add_argument(
        "--writer-queue-size",
        type=int,
        default=2_000,
        help="Max queue items between embedder and writer (default: 2000)",
    )
    parser.add_argument(
        "--writer-batch-size",
        type=int,
        default=1_000,
        help="Embeddings to flush to Milvus per batch (default: 1000)",
    )
    parser.add_argument(
        "--status-batch-size",
        type=int,
        default=1_000,
        help="Embeddings to trigger status flush (default: 1000)",
    )
    parser.add_argument(
        "--huggingface-token",
        type=str,
        default=None,
        help="Optional HuggingFace token (otherwise uses cached login)",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable Rich progress bar output",
    )

    args = parser.parse_args(argv)

    try:
        devices = parse_devices(args.devices)
    except ValueError as exc:
        parser.error(str(exc))

    driver = get_async_driver()
    milvus_client = get_async_milvus_client()

    progress_columns = [
        SpinnerColumn(),
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
    ]

    try:
        if args.no_progress:
            await run_embedding_job(
                neo4j_driver=driver,
                milvus_client=milvus_client,
                collection_name=args.collection,
                devices=devices,
                dtype=args.dtype,
                max_padded_tokens=args.max_padded_tokens,
                max_length=args.max_length,
                source_chunk_size=args.source_chunk_size,
                read_batch_size=args.read_batch_size,
                writer_queue_size=args.writer_queue_size,
                writer_batch_size=args.writer_batch_size,
                status_batch_size=args.status_batch_size,
                huggingface_token=args.huggingface_token,
                progress=None,
                progress_task_id=None,
            )
        else:
            with Progress(*progress_columns) as progress:
                task_id = progress.add_task("embedding", total=None)
                await run_embedding_job(
                    neo4j_driver=driver,
                    milvus_client=milvus_client,
                    collection_name=args.collection,
                    devices=devices,
                    dtype=args.dtype,
                    max_padded_tokens=args.max_padded_tokens,
                    max_length=args.max_length,
                    source_chunk_size=args.source_chunk_size,
                    read_batch_size=args.read_batch_size,
                    writer_queue_size=args.writer_queue_size,
                    writer_batch_size=args.writer_batch_size,
                    status_batch_size=args.status_batch_size,
                    huggingface_token=args.huggingface_token,
                    progress=progress,
                    progress_task_id=task_id,
                )
    finally:
        await milvus_client.close()
        await driver.close()


if __name__ == "__main__":
    asyncio.run(main())
