"""Pipeline for embedding pending proteins from Neo4j and storing in Milvus.

Streaming pipeline with length-sorted batching:
    [Neo4j bulk fetch] → [Sort by length] → [ESM2 Embedder] → [Milvus + Status Update]

Memory bounded via prefetch size and embedder backpressure.
Sequences are sorted by length to batch similar-sized sequences together,
avoiding ESM2 rotary embedding cache mismatches.

Example:
    from pyeed.db.neo4j import get_async_driver
    from pyeed.db.milvus import get_async_milvus_client
    from pyeed.embed.embed_pending_proteins import embed_pending_proteins

    async def main():
        driver = get_async_driver()
        client = get_async_milvus_client()

        stats = await embed_pending_proteins(
            driver=driver,
            client=client,
            collection_name="proteins",
        )
        print(stats)  # {'completed': 9950, 'failed': 50}
"""

from __future__ import annotations

import asyncio
from collections.abc import AsyncIterator
from typing import Literal

from loguru import logger
from neo4j import AsyncDriver
from pymilvus import AsyncMilvusClient
from rich.progress import (
    MofNCompleteColumn,
    Progress,
    ProgressColumn,
    SpinnerColumn,
    Task,
    Text,
    TextColumn,
    TimeElapsedColumn,
)

from pyeed.db.milvus import initialize_collection_from_dict, insert
from pyeed.embed.esm2 import ESM2Embedder

type DType = Literal["float16", "float32"]


class IterPerSecColumn(ProgressColumn):
    """Show iterations per second."""

    def render(self, task: Task) -> Text:
        if task.speed is None:
            return Text("- it/s")
        return Text(f"{task.speed:.0f} it/s")


# =============================================================================
# Database Operations
# =============================================================================


async def update_protein_status(
    driver: AsyncDriver,
    protein_ids: list[str],
    status: Literal["pending", "in_progress", "complete", "failed"],
) -> None:
    """Update embedding_status for proteins using direct Cypher.

    Uses a simple SET operation that only modifies the embedding_status field,
    preserving all other protein data.

    Args:
        driver: Neo4j async driver.
        protein_ids: List of protein IDs to update.
        status: New embedding status value.
    """
    if not protein_ids:
        return

    query = """
    UNWIND $ids AS id
    MATCH (p:Protein {id: id})
    SET p.embedding_status = $status
    """
    async with driver.session() as session:
        await session.run(query, ids=protein_ids, status=status)

    logger.debug(f"Updated {len(protein_ids)} proteins to status '{status}'")


async def fetch_pending_batch(
    driver: AsyncDriver,
    limit: int,
    max_seq_length: int | None = None,
) -> list[tuple[str, str]]:
    conditions = ["p.embedding_status = 'pending'"]
    if max_seq_length:
        conditions.append("p.seq_length <= $max_seq_length")

    where = " AND ".join(conditions)

    query = """
    CALL () {
    MATCH (p:Protein)
    WHERE p.embedding_status = 'pending'
    RETURN p.id AS id, p.sequence AS sequence
    LIMIT $limit
    }
    WITH id, sequence
    MATCH (p:Protein {id: id})
    SET p.embedding_status = 'in_progress'
    RETURN id, sequence;
    """

    params: dict = {"limit": limit}
    if max_seq_length:
        params["max_seq_length"] = max_seq_length

    async with driver.session() as session:
        result = await session.run(query, params)
        records = await result.data()

    return [(r["id"], r["sequence"]) for r in records]


async def stream_pending_proteins(
    driver: AsyncDriver,
    batch_size: int = 64,
    prefetch_size: int = 1024,
    max_seq_length: int | None = None,
) -> AsyncIterator[tuple[list[str], list[str]]]:
    """Stream pending proteins as length-sorted (ids, sequences) batches.

    Fetches proteins in large batches, sorts by sequence length locally,
    then yields smaller batches of similar-length sequences.

    Args:
        driver: Neo4j async driver.
        batch_size: Number of proteins per yielded batch.
        prefetch_size: Number of proteins to fetch and sort at once.
        max_seq_length: Optional maximum sequence length filter.

    Yields:
        Tuples of (ids_list, sequences_list) sorted by sequence length.
    """
    prefetch_num = 0

    while True:
        # Fetch large batch from Neo4j
        logger.debug(f"Fetching {prefetch_size} pending proteins")
        proteins = await fetch_pending_batch(driver, prefetch_size, max_seq_length)
        logger.debug(f"Fetched {len(proteins)} proteins")

        if not proteins:
            logger.info("No more pending proteins")
            break

        prefetch_num += 1

        # Sort by sequence length to group similar-sized sequences
        proteins_sorted = sorted(proteins, key=lambda p: len(p[1]))

        min_len = len(proteins_sorted[0][1])
        max_len = len(proteins_sorted[-1][1])
        logger.info(
            f"Prefetch {prefetch_num}: {len(proteins)} proteins (length {min_len}-{max_len})"
        )

        # Yield small batches of similar-length sequences
        for i in range(0, len(proteins_sorted), batch_size):
            batch = proteins_sorted[i : i + batch_size]
            ids = [p[0] for p in batch]
            sequences = [p[1] for p in batch]

            batch_min = len(sequences[0])
            batch_max = len(sequences[-1])
            logger.debug(f"Batch: {len(ids)} proteins (length {batch_min}-{batch_max})")

            yield ids, sequences

            # Allow event loop to process other tasks
            await asyncio.sleep(0)


async def flush_accumulated_batch(
    driver: AsyncDriver,
    client: AsyncMilvusClient,
    collection_name: str,
    accumulated_ids: list[str],
    accumulated_records: list[dict],
) -> dict[str, int]:
    """Flush accumulated records to Milvus and update status.

    Returns:
        {'completed': n_completed, 'failed': n_failed}
    """
    if not accumulated_records:
        return {"completed": 0, "failed": 0}

    completed = 0
    failed = 0

    try:
        await insert(
            collection_name=collection_name,
            records=accumulated_records,
            client=client,
        )
        await update_protein_status(driver, accumulated_ids, "complete")
        completed = len(accumulated_ids)

        logger.info(f"Flushed {len(accumulated_ids)} embeddings to Milvus")

    except Exception as e:
        logger.error(f"Insert failed for {len(accumulated_ids)} proteins: {e}")
        await update_protein_status(driver, accumulated_ids, "failed")
        failed = len(accumulated_ids)

    return {"completed": completed, "failed": failed}


# =============================================================================
# Main Pipeline
# =============================================================================


async def embed_pending_proteins(
    driver: AsyncDriver,
    client: AsyncMilvusClient,
    collection_name: str,
    batch_size: int = 64,
    prefetch_size: int = 12800,
    insert_batch_size: int = 1000,
    model_name: str = "facebook/esm2_t33_650M_UR50D",
    dtype: DType = "float16",
    max_length: int = 1024,
    devices: list[int] | None = None,
    max_seq_length: int | None = None,
    max_flush_in_flight: int = 4,
) -> dict[str, int]:
    if devices is None:
        devices = [0]

    stats = {"completed": 0, "failed": 0}
    collection_initialized = collection_name in await client.list_collections()

    accumulated_ids: list[str] = []
    accumulated_records: list[dict] = []
    flush_tasks: set[asyncio.Task[dict[str, int]]] = set()

    logger.info(
        f"Starting embedding pipeline (batch_size={batch_size}, "
        f"prefetch_size={prefetch_size}, insert_batch_size={insert_batch_size}, "
        f"devices={devices})"
    )

    with Progress(
        SpinnerColumn(),
        TextColumn("[bold blue]{task.description}"),
        IterPerSecColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("[cyan]Embedding proteins", total=None)

        try:
            async with ESM2Embedder(
                model_name=model_name,
                dtype=dtype,
                max_length=max_length,
                devices=devices,
            ) as embedder:
                protein_stream = stream_pending_proteins(
                    driver, batch_size, prefetch_size, max_seq_length
                )

                pending_batch_ids: list[list[str]] = []

                async def tracked_stream():
                    async for ids, sequences in protein_stream:
                        pending_batch_ids.append(ids)
                        yield ids, sequences

                async for batch_results in embedder.embed_stream(tracked_stream()):
                    original_ids = pending_batch_ids.pop(0) if pending_batch_ids else []
                    if not original_ids:
                        continue

                    embedded_ids = {pid for pid, _ in batch_results} if batch_results else set()
                    skipped_ids = [pid for pid in original_ids if pid not in embedded_ids]

                    if skipped_ids:
                        logger.warning(
                            f"Marking {len(skipped_ids)} skipped proteins as failed "
                            "(exceeded max_length)"
                        )
                        await update_protein_status(driver, skipped_ids, "failed")
                        stats["failed"] += len(skipped_ids)
                        progress.update(task, advance=len(skipped_ids))

                    if not batch_results:
                        continue

                    ids = [pid for pid, _ in batch_results]
                    records = [{"id": pid, "embedding": emb} for pid, emb in batch_results]

                    if not collection_initialized:
                        await initialize_collection_from_dict(
                            collection_name=collection_name,
                            primary_field_name="id",
                            record=records[0],
                            client=client,
                        )
                        collection_initialized = True
                        logger.info(f"Initialized collection '{collection_name}'")

                    accumulated_ids.extend(ids)
                    accumulated_records.extend(records)

                    logger.debug(
                        f"Accumulated {len(accumulated_records)}/{insert_batch_size} records"
                    )
                    # Progress is about proteins processed, not flushes:
                    progress.update(task, advance=len(records))

                    if len(accumulated_records) >= insert_batch_size:
                        # snapshot current batch for flushing
                        ids_to_flush = accumulated_ids
                        records_to_flush = accumulated_records
                        accumulated_ids = []
                        accumulated_records = []

                        flush_task = asyncio.create_task(
                            flush_accumulated_batch(
                                driver,
                                client,
                                collection_name,
                                ids_to_flush,
                                records_to_flush,
                            )
                        )
                        flush_tasks.add(flush_task)

                        # limit number of concurrent flushes
                        if len(flush_tasks) >= max_flush_in_flight:
                            logger.debug(f"Waiting for {len(flush_tasks)} flush tasks to complete")
                            done, flush_tasks = await asyncio.wait(
                                flush_tasks, return_when=asyncio.FIRST_COMPLETED
                            )
                            for t in done:
                                delta = t.result()
                                stats["completed"] += delta["completed"]
                                stats["failed"] += delta["failed"]

        except (KeyboardInterrupt, asyncio.CancelledError) as e:
            logger.warning(f"Pipeline interrupted: {type(e).__name__}")
        finally:
            # flush remaining accumulated records
            if accumulated_records:
                logger.info(f"Scheduling final flush of {len(accumulated_records)} records")
                flush_task = asyncio.create_task(
                    flush_accumulated_batch(
                        driver,
                        client,
                        collection_name,
                        accumulated_ids,
                        accumulated_records,
                    )
                )
                flush_tasks.add(flush_task)

            # wait for all outstanding flushes and merge stats
            if flush_tasks:
                done, _ = await asyncio.wait(flush_tasks)
                for t in done:
                    delta = t.result()
                    stats["completed"] += delta["completed"]
                    stats["failed"] += delta["failed"]

    logger.info(f"Pipeline complete: {stats['completed']} completed, {stats['failed']} failed")
    return stats


# =============================================================================
# Entry Point
# =============================================================================


if __name__ == "__main__":
    import argparse
    import sys

    from rich import print as rprint

    from pyeed.db.milvus import get_async_milvus_client
    from pyeed.db.neo4j import get_async_driver

    def parse_devices(devices_str: str) -> list[int]:
        """Parse comma-separated device IDs into list of integers.

        Args:
            devices_str: Comma-separated device IDs, e.g., "0,1,2" or "2,0".

        Returns:
            List of device IDs as integers.

        Raises:
            ValueError: If device IDs cannot be parsed as integers.
        """
        if not devices_str:
            return [0]
        try:
            return [int(d.strip()) for d in devices_str.split(",") if d.strip()]
        except ValueError as e:
            msg = (
                f"Invalid device IDs format: {devices_str}. "
                "Use comma-separated integers, e.g., '0,1,2'"
            )
            raise ValueError(msg) from e

    async def main() -> None:
        parser = argparse.ArgumentParser(
            description="Embed pending proteins from Neo4j and store in Milvus",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Examples:
  # Use default settings (device 0, collection 'proteins')
  python -m pyeed.embed.embed_pending_proteins

  # Use multiple GPUs (devices 2 and 0)
  python -m pyeed.embed.embed_pending_proteins --devices 2,0

  # Custom collection and batch sizes
  python -m pyeed.embed.embed_pending_proteins --collection my_proteins --batch-size 128

  # Limit to sequences <= 512 amino acids
  python -m pyeed.embed.embed_pending_proteins --max-seq-length 512
            """,
        )

        parser.add_argument(
            "--collection",
            type=str,
            default="proteins",
            help="Name of Milvus collection to store embeddings (default: 'proteins')",
        )
        parser.add_argument(
            "--devices",
            type=str,
            default="0",
            help="Comma-separated CUDA device IDs, e.g., '0,1,2' or '2,0' (default: '0')",
        )
        parser.add_argument(
            "--batch-size",
            type=int,
            default=64,
            help="Number of proteins per embedding batch (default: 64)",
        )
        parser.add_argument(
            "--prefetch-size",
            type=int,
            default=12800,
            help="Number of proteins to fetch and sort at once (default: 12800)",
        )
        parser.add_argument(
            "--insert-batch-size",
            type=int,
            default=1024,
            help="Number of embeddings to accumulate before inserting (default: 1024)",
        )
        parser.add_argument(
            "--model-name",
            type=str,
            default="facebook/esm2_t33_650M_UR50D",
            help="HuggingFace model identifier for ESM2 (default: 'facebook/esm2_t33_650M_UR50D')",
        )
        parser.add_argument(
            "--dtype",
            type=str,
            choices=["float16", "float32"],
            default="float16",
            help="Data type for model computation and embeddings (default: 'float16')",
        )
        parser.add_argument(
            "--max-length",
            type=int,
            default=1024,
            help="Maximum token length for ESM2 (default: 1024)",
        )
        parser.add_argument(
            "--max-seq-length",
            type=int,
            default=None,
            help="Optional maximum amino acid sequence length filter (default: None)",
        )
        parser.add_argument(
            "--max-flush-in-flight",
            type=int,
            default=4,
            help="Maximum number of concurrent flush operations (default: 4)",
        )

        args = parser.parse_args()

        try:
            devices = parse_devices(args.devices)
        except ValueError as e:
            rprint(f"[red]Error:[/red] {e}")
            sys.exit(1)

        driver = get_async_driver()
        client = get_async_milvus_client()

        rprint(
            f"[cyan]Starting embedding pipeline...[/cyan]\n"
            f"  Collection: {args.collection}\n"
            f"  Devices: {devices}\n"
            f"  Batch size: {args.batch_size}\n"
            f"  Prefetch size: {args.prefetch_size}\n"
            f"  Insert batch size: {args.insert_batch_size}\n"
            f"  Model: {args.model_name}\n"
            f"  Dtype: {args.dtype}\n"
            f"  Max length: {args.max_length}"
        )
        if args.max_seq_length:
            rprint(f"  Max sequence length: {args.max_seq_length}")

        try:
            stats = await embed_pending_proteins(
                driver=driver,
                client=client,
                collection_name=args.collection,
                batch_size=args.batch_size,
                prefetch_size=args.prefetch_size,
                insert_batch_size=args.insert_batch_size,
                model_name=args.model_name,
                dtype=args.dtype,
                max_length=args.max_length,
                devices=devices,
                max_seq_length=args.max_seq_length,
                max_flush_in_flight=args.max_flush_in_flight,
            )
            rprint("[green]Done![/green]")
            rprint(stats)
        finally:
            await client.close()
            await driver.close()

    asyncio.run(main())
