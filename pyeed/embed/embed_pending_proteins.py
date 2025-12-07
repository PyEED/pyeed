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
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

from pyeed.db.milvus import initialize_collection_from_dict, insert
from pyeed.embed.esm2 import ESM2Embedder

type DType = Literal["float16", "float32"]


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
    """Fetch a batch of pending proteins from Neo4j.

    Args:
        driver: Neo4j async driver.
        limit: Maximum number of proteins to fetch.
        max_seq_length: Optional maximum sequence length filter.

    Returns:
        List of (id, sequence) tuples.
    """
    conditions = ["p.sequence IS NOT NULL"]
    if max_seq_length:
        conditions.append("p.seq_length <= $max_seq_length")

    where = " AND ".join(conditions)
    query = f"""
    MATCH (p:Protein)
    WHERE {where}
      AND (p.embedding_status IS NULL OR p.embedding_status = 'pending')
    RETURN p.id AS id, p.sequence AS sequence
    LIMIT $limit
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
    prefetch_size: int = 12800,
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
        proteins = await fetch_pending_batch(driver, prefetch_size, max_seq_length)

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

            # Mark as in_progress BEFORE yielding
            await update_protein_status(driver, ids, "in_progress")

            yield ids, sequences

            # Allow event loop to process other tasks
            await asyncio.sleep(0)


async def flush_accumulated_batch(
    driver: AsyncDriver,
    client: AsyncMilvusClient,
    collection_name: str,
    accumulated_ids: list[str],
    accumulated_records: list[dict],
    stats: dict[str, int],
) -> None:
    """Flush accumulated records to Milvus and update status.

    Args:
        driver: Neo4j async driver.
        client: AsyncMilvusClient instance.
        collection_name: Name of Milvus collection.
        accumulated_ids: List of protein IDs in batch.
        accumulated_records: List of records to insert.
        stats: Statistics dict to update.
    """
    if not accumulated_records:
        return

    try:
        await insert(
            collection_name=collection_name,
            records=accumulated_records,
            client=client,
        )
        await update_protein_status(driver, accumulated_ids, "complete")
        stats["completed"] += len(accumulated_ids)

        logger.info(
            f"Flushed {len(accumulated_ids)} embeddings to Milvus (total: {stats['completed']})"
        )

    except Exception as e:
        logger.error(f"Insert failed for {len(accumulated_ids)} proteins: {e}")
        await update_protein_status(driver, accumulated_ids, "failed")
        stats["failed"] += len(accumulated_ids)


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
    device: int = 0,
    max_seq_length: int | None = None,
) -> dict[str, int]:
    """Embed pending proteins: stream → embed → insert → update status.

    Fetches proteins in large batches, sorts by length to avoid ESM2 cache
    mismatches, then streams through the embedder with built-in backpressure.
    Accumulates embeddings and inserts in batches to reduce database overhead.

    Args:
        driver: Neo4j async driver.
        client: AsyncMilvusClient instance.
        collection_name: Name of Milvus collection to store embeddings.
        batch_size: Number of proteins per embedding batch.
        prefetch_size: Number of proteins to fetch and sort at once.
        insert_batch_size: Number of embeddings to accumulate before inserting.
        model_name: HuggingFace model identifier for ESM2.
        dtype: Data type for model computation and embeddings.
        max_length: Maximum token length for ESM2.
        device: GPU device index (0, 1, ...) or -1 for all GPUs.
        max_seq_length: Optional maximum amino acid sequence length filter.

    Returns:
        Statistics dict with 'completed' and 'failed' counts.
    """
    stats = {"completed": 0, "failed": 0}
    collection_initialized = collection_name in await client.list_collections()

    # Accumulators for batched inserts
    accumulated_ids: list[str] = []
    accumulated_records: list[dict] = []

    logger.info(
        f"Starting embedding pipeline (batch_size={batch_size}, "
        f"prefetch_size={prefetch_size}, insert_batch_size={insert_batch_size}, "
        f"device={device})"
    )

    with Progress(
        SpinnerColumn(),
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("[cyan]Embedding proteins", total=None)

        try:
            async with ESM2Embedder(
                model_name=model_name,
                dtype=dtype,
                max_length=max_length,
                device=device,
            ) as embedder:
                # Stream length-sorted proteins into embedder
                protein_stream = stream_pending_proteins(
                    driver, batch_size, prefetch_size, max_seq_length
                )

                # Track IDs as they enter the embedding pipeline
                pending_batch_ids: list[list[str]] = []

                async def tracked_stream() -> AsyncIterator[tuple[list[str], list[str]]]:
                    """Stream that tracks original IDs."""
                    async for ids, sequences in protein_stream:
                        pending_batch_ids.append(ids)
                        yield ids, sequences

                async for batch_results in embedder.embed_stream(tracked_stream()):
                    # Get the original IDs for this batch
                    original_ids = pending_batch_ids.pop(0) if pending_batch_ids else []

                    if not original_ids:
                        continue

                    # Extract successfully embedded IDs
                    embedded_ids = {pid for pid, _ in batch_results} if batch_results else set()

                    # Identify skipped IDs (marked in_progress but not embedded)
                    skipped_ids = [pid for pid in original_ids if pid not in embedded_ids]

                    if skipped_ids:
                        logger.warning(
                            f"Marking {len(skipped_ids)} skipped proteins as failed "
                            f"(exceeded max_length)"
                        )
                        await update_protein_status(driver, skipped_ids, "failed")
                        stats["failed"] += len(skipped_ids)
                        progress.update(task, advance=len(skipped_ids))

                    if not batch_results:
                        continue

                    ids = [pid for pid, _ in batch_results]
                    records = [{"id": pid, "embedding": emb} for pid, emb in batch_results]

                    # Initialize collection lazily on first batch
                    if not collection_initialized:
                        await initialize_collection_from_dict(
                            collection_name=collection_name,
                            primary_field_name="id",
                            record=records[0],
                            client=client,
                        )
                        collection_initialized = True
                        logger.info(f"Initialized collection '{collection_name}'")

                    # Accumulate records
                    accumulated_ids.extend(ids)
                    accumulated_records.extend(records)

                    logger.debug(
                        f"Accumulated {len(accumulated_records)}/{insert_batch_size} records"
                    )
                    progress.update(task, advance=len(records))

                    # Flush when we reach the insert batch size
                    if len(accumulated_records) >= insert_batch_size:
                        await flush_accumulated_batch(
                            driver,
                            client,
                            collection_name,
                            accumulated_ids,
                            accumulated_records,
                            stats,
                        )

                        # Clear accumulators
                        accumulated_ids = []
                        accumulated_records = []

        except (KeyboardInterrupt, asyncio.CancelledError) as e:
            logger.warning(f"Pipeline interrupted: {type(e).__name__}")
        finally:
            # Always flush remaining records on exit (normal or interrupted)
            if accumulated_records:
                logger.info(f"Flushing final batch of {len(accumulated_records)} records")
                await flush_accumulated_batch(
                    driver,
                    client,
                    collection_name,
                    accumulated_ids,
                    accumulated_records,
                    stats,
                )
                progress.update(task, advance=len(accumulated_ids))

    logger.info(f"Pipeline complete: {stats['completed']} completed, {stats['failed']} failed")
    return stats


# =============================================================================
# Entry Point
# =============================================================================


if __name__ == "__main__":
    from rich import print as rprint

    from pyeed.db.milvus import get_async_milvus_client
    from pyeed.db.neo4j import get_async_driver

    async def main() -> None:
        driver = get_async_driver()
        client = get_async_milvus_client()

        try:
            stats = await embed_pending_proteins(
                driver=driver,
                client=client,
                collection_name="proteins",
                batch_size=64,
                prefetch_size=12800,
                insert_batch_size=1000,
                dtype="float16",
                device=0,
            )
            rprint(f"[green]Done:[/green] {stats}")
        finally:
            await client.close()
            await driver.close()

    asyncio.run(main())
