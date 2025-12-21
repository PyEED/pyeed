from __future__ import annotations

import argparse
import asyncio
from collections.abc import AsyncIterator
from time import monotonic
from typing import Final, Literal

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
from .embedder import ESM2Processor
from .types import EmbeddingRecord, SequenceItem, TorchDTypeName

__all__ = [
    "run_embedding_job",
    "stream_pending_proteins",
    "update_protein_status",
]

type StatusLiteral = Literal["pending", "in_progress", "complete", "failed"]

_PENDING_QUERY: Final[str] = """
MATCH (p:Protein {embedding_status: 'pending'})
RETURN p.id AS id, p.sequence AS sequence
"""


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
) -> AsyncIterator[list[SequenceItem]]:
    """Stream pending protein ids and sequences in bounded batches.

    Keeps the Neo4j session open for streaming to avoid loading all rows into
    memory while still yielding manageable batches for the embedder.
    Marks yielded proteins as ``in_progress`` when mark_in_progress is True.
    """
    async with driver.session() as session:
        result = await session.run(_PENDING_QUERY)

        buffer: list[SequenceItem] = []
        async for record in result:
            protein_id = record["id"]
            sequence = record["sequence"]
            buffer.append(SequenceItem(id=protein_id, sequence=sequence))

            if len(buffer) >= batch_size:
                logger.debug(
                    "Yielding pending batch",
                    extra={"batch_size": len(buffer)},
                )
                if mark_in_progress:
                    await update_protein_status(driver, [item.id for item in buffer], "in_progress")
                yield buffer
                buffer = []

        if buffer:
            logger.debug(
                "Yielding final pending batch",
                extra={"batch_size": len(buffer)},
            )
            if mark_in_progress:
                await update_protein_status(driver, [item.id for item in buffer], "in_progress")
            yield buffer


async def run_embedding_job(
    *,
    neo4j_driver: AsyncDriver,
    milvus_client: AsyncMilvusClient,
    collection_name: str,
    devices: list[int] | None = None,
    dtype: TorchDTypeName = "float16",
    max_padded_tokens: int = 30_000,
    max_length: int = 1_024,
    source_chunk_size: int = 2_048,
    read_batch_size: int = 5_000,
    writer_queue_size: int = 10_000,
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
    logger.info("Embedding job starting")
    if progress is not None and progress_task_id is not None:
        progress.update(progress_task_id, total=None)

    out_q: asyncio.Queue[EmbeddingRecord | None] = asyncio.Queue(writer_queue_size)
    collection_initialized = collection_name in await milvus_client.list_collections()

    async def flush(batch: list[EmbeddingRecord]) -> None:
        nonlocal collection_initialized
        if not batch:
            return
        batch_count = len(batch)
        batch_ids = [record.id for record in batch]
        start = monotonic()
        records = [{"id": record.id, "embedding": record.vector} for record in batch]

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

            insert_start = monotonic()
            await insert(
                collection_name=collection_name,
                records=records,
                client=milvus_client,
            )
            insert_elapsed = monotonic() - insert_start

            status_start = monotonic()
            await update_protein_status(neo4j_driver, batch_ids, "complete")
            status_elapsed = monotonic() - status_start

            total_elapsed = round(monotonic() - start, 3)
            insert_elapsed_rounded = round(insert_elapsed, 3)
            status_elapsed_rounded = round(status_elapsed, 3)
            logger.debug(
                f"Writer flushed batch in {total_elapsed}s",
                extra={
                    "items": batch_count,
                    "total_elapsed_s": total_elapsed,
                    "insert_elapsed_s": insert_elapsed_rounded,
                    "status_elapsed_s": status_elapsed_rounded,
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

    async def flush_with_semaphore(
        batch: list[EmbeddingRecord],
        sem: asyncio.Semaphore,
    ) -> None:
        async with sem:
            await flush(batch)  # your existing flush logic

    async def writer() -> None:
        # Allow up to 3 concurrent flushes
        sem = asyncio.Semaphore(3)
        buffer: list[EmbeddingRecord] = []
        flush_tasks: list[asyncio.Task] = []

        while True:
            item = await out_q.get()
            try:
                if item is None:
                    break
                buffer.append(item)

                if len(buffer) >= writer_batch_size:
                    # Start flush in background
                    task = asyncio.create_task(flush_with_semaphore(buffer.copy(), sem))
                    flush_tasks.append(task)
                    buffer.clear()
            finally:
                out_q.task_done()

        # Flush remaining and wait for all tasks
        if buffer:
            await flush(buffer)
        await asyncio.gather(*flush_tasks)

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
        async for record in processor.stream_work(chunks):
            put_start = monotonic()
            await out_q.put(record)
            waited = monotonic() - put_start
            if waited > 0.1:
                logger.debug(
                    "Backpressure on embedding queue",
                    extra={
                        "waited_s": round(waited, 3),
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
        default=50_000,
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
        default=2000,
        help="Max queue items between embedder and writer (default: 2000)",
    )
    parser.add_argument(
        "--writer-batch-size",
        type=int,
        default=400,
        help="Embeddings to flush to Milvus per batch (default: 2000)",
    )
    parser.add_argument(
        "--status-batch-size",
        type=int,
        default=400,
        help="Embeddings to trigger status flush (default: 2000)",
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
