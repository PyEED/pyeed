"""Pipeline for embedding proteins from Neo4j and storing in Milvus.

This module provides functions to:
1. Fetch proteins with sequences from Neo4j in batches
2. Filter by sequence length (optional)
3. Embed using ESM2Embedder with backpressure
4. Store embeddings in Milvus

Example:
    from pyeed.db.neo4j import get_async_driver
    from pyeed.db.milvus import get_async_milvus_client
    from pyeed.embed.embed_proteins import embed_proteins, get_embeddable_proteins

    async def main():
        driver = get_async_driver()
        client = get_async_milvus_client()

        # Get proteins not yet embedded
        protein_ids = await get_embeddable_proteins(
            driver, client, "proteins", "id"
        )

        # Embed and store
        stats = await embed_proteins(
            driver=driver,
            client=client,
            collection_name="proteins",
            ids=protein_ids,
            dtype="float16",
            device=0,
        )
        print(stats)
"""

from __future__ import annotations

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
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from pyeed.db.milvus import initialize_collection_from_dict, insert
from pyeed.embed.esm2 import ESM2Embedder

type DType = Literal["float16", "float32"]


def create_progress() -> Progress:
    """Create a Rich progress bar for protein embedding."""
    return Progress(
        SpinnerColumn(),
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        refresh_per_second=4,
    )


async def fetch_sequence_batches(
    driver: AsyncDriver,
    batch_size: int = 1000,
    ids: list[str] | None = None,
    max_seq_length: int | None = None,
) -> AsyncIterator[list[tuple[str, str]]]:
    """Fetch proteins with sequences from Neo4j in batches.

    Queries proteins that have a non-null sequence attribute. Optionally filters
    by a specific list of protein IDs and/or maximum sequence length.

    Args:
        driver: Neo4j async driver.
        batch_size: Number of proteins to fetch per query.
        ids: Optional list of protein IDs to filter by. If None, fetches all
            embeddable proteins.
        max_seq_length: Optional maximum sequence length filter. If None, no
            length filtering is applied.

    Yields:
        Batches of (id, sequence) tuples for valid proteins.
    """
    # Build query based on filters
    conditions = ["p.sequence IS NOT NULL"]
    if ids:
        conditions.append("p.id IN $ids")
    if max_seq_length:
        conditions.append("p.seq_length <= $max_seq_length")

    where_clause = " AND ".join(conditions)
    query = f"""
    MATCH (p:Protein)
    WHERE {where_clause}
    RETURN p.id AS id, p.sequence AS sequence
    SKIP $offset LIMIT $limit
    """

    offset = 0
    while True:
        async with driver.session() as session:
            params: dict = {"offset": offset, "limit": batch_size}
            if ids:
                params["ids"] = ids
            if max_seq_length:
                params["max_seq_length"] = max_seq_length

            result = await session.run(query, params)
            records = await result.data()

        if not records:
            break

        batch = [(r["id"], r["sequence"]) for r in records]
        logger.debug(f"Fetched {len(batch)} proteins from Neo4j", extra={"offset": offset})
        yield batch

        if len(records) < batch_size:
            break

        offset += batch_size


async def rebatch_for_embedder(
    source: AsyncIterator[list[tuple[str, str]]],
    batch_size: int = 32,
    fetched_counter: list[int] | None = None,
) -> AsyncIterator[tuple[list[str], list[str]]]:
    """Rebatch proteins into fixed-size batches for the embedder.

    Takes variable-sized batches from source and yields fixed-size
    (ids, sequences) tuples suitable for the embedder.

    Args:
        source: Async iterator yielding batches of (id, sequence) tuples.
        batch_size: Target batch size for embedder.
        fetched_counter: Optional mutable list to track total fetched count.

    Yields:
        Tuples of (ids_list, sequences_list) with up to batch_size items.
    """
    buffer_ids: list[str] = []
    buffer_sequences: list[str] = []

    async for batch in source:
        for protein_id, sequence in batch:
            buffer_ids.append(protein_id)
            buffer_sequences.append(sequence)
            if fetched_counter is not None:
                fetched_counter[0] += 1

            if len(buffer_ids) >= batch_size:
                yield buffer_ids[:batch_size], buffer_sequences[:batch_size]
                buffer_ids = buffer_ids[batch_size:]
                buffer_sequences = buffer_sequences[batch_size:]

    # Yield remainder
    if buffer_ids:
        yield buffer_ids, buffer_sequences


async def fetch_all_embeddable_protein_ids(
    driver: AsyncDriver,
    max_seq_length: int | None = None,
) -> list[str]:
    """Fetch all protein IDs with valid sequences from Neo4j.

    Args:
        driver: Neo4j async driver.
        max_seq_length: Optional maximum sequence length filter.

    Returns:
        List of protein IDs with non-null sequences.
    """
    conditions = ["p.sequence IS NOT NULL"]
    if max_seq_length:
        conditions.append("p.seq_length <= $max_seq_length")

    where_clause = " AND ".join(conditions)
    query = f"""
    MATCH (p:Protein)
    WHERE {where_clause}
    RETURN p.id AS id
    """

    ids: list[str] = []
    async with driver.session() as session:
        params = {"max_seq_length": max_seq_length} if max_seq_length else {}
        result = await session.run(query, params)
        async for record in result:
            if record["id"] is not None:
                ids.append(record["id"])

    return ids


async def count_embeddable_proteins(
    driver: AsyncDriver,
    ids: list[str] | None = None,
    max_seq_length: int | None = None,
) -> int:
    """Count proteins with valid sequences for embedding.

    Args:
        driver: Neo4j async driver.
        ids: Optional list of protein IDs to count.
        max_seq_length: Optional maximum sequence length filter.

    Returns:
        Total count of proteins with non-null sequences.
    """
    conditions = ["p.sequence IS NOT NULL"]
    if ids:
        conditions.append("p.id IN $ids")
    if max_seq_length:
        conditions.append("p.seq_length <= $max_seq_length")

    where_clause = " AND ".join(conditions)
    query = f"""
    MATCH (p:Protein)
    WHERE {where_clause}
    RETURN count(p) AS count
    """

    async with driver.session() as session:
        params: dict = {}
        if ids:
            params["ids"] = ids
        if max_seq_length:
            params["max_seq_length"] = max_seq_length

        result = await session.run(query, params)
        record = await result.single()
        return record["count"] if record else 0


async def get_embeddable_proteins(
    driver: AsyncDriver,
    client: AsyncMilvusClient,
    collection_name: str,
    primary_field_name: str,
    max_seq_length: int | None = None,
) -> list[str]:
    """Return protein IDs that need embedding (not yet in Milvus).

    Fetches all proteins with valid sequences from Neo4j and filters out
    those that already exist in the Milvus collection.

    Args:
        driver: Neo4j async driver.
        client: Milvus async client.
        collection_name: Name of the Milvus collection to check.
        primary_field_name: Name of the primary key field in Milvus.
        max_seq_length: Optional maximum sequence length filter.

    Returns:
        List of protein IDs that need embedding (not in collection).
        Returns all embeddable protein IDs if collection doesn't exist.
    """
    # Fetch all embeddable protein IDs from Neo4j
    all_ids = await fetch_all_embeddable_protein_ids(driver, max_seq_length)

    if not all_ids:
        return []

    # Check if collection exists
    existing_collections = await client.list_collections()
    if collection_name not in existing_collections:
        logger.info(
            f"Collection '{collection_name}' does not exist, returning all {len(all_ids)} protein IDs"
        )
        return all_ids

    # Get existing IDs from Milvus
    try:
        existing_records = await client.get(
            collection_name=collection_name,
            ids=all_ids,
            output_fields=[primary_field_name],
        )
        existing_ids = {record[primary_field_name] for record in existing_records}
    except Exception as e:
        logger.warning(
            f"Failed to check existing IDs in Milvus: {e}",
            extra={"collection": collection_name},
        )
        return all_ids

    logger.info(
        f"Found {len(all_ids)} embeddable proteins in Neo4j, "
        f"{len(existing_ids)} already in Milvus, "
        f"{len(all_ids) - len(existing_ids)} need embedding",
        extra={
            "neo4j_count": len(all_ids),
            "milvus_count": len(existing_ids),
            "to_embed": len(all_ids) - len(existing_ids),
        },
    )

    # Return IDs that need embedding
    return [protein_id for protein_id in all_ids if protein_id not in existing_ids]


async def embed_proteins(
    driver: AsyncDriver,
    client: AsyncMilvusClient,
    collection_name: str,
    graph_db_batch_size: int = 1000,
    embed_batch_size: int = 32,
    model_name: str = "facebook/esm2_t33_650M_UR50D",
    dtype: DType = "float16",
    max_length: int = 1024,
    device: int = -1,
    show_progress: bool = True,
    ids: list[str] | None = None,
    max_seq_length: int | None = None,
) -> dict[str, int]:
    """Embed proteins from Neo4j and store embeddings in Milvus.

    Fetches proteins with sequences from Neo4j, embeds them using the
    ESM2Embedder, and inserts the embeddings into Milvus. Creates
    the collection automatically if it doesn't exist.

    Args:
        driver: Neo4j async driver.
        client: Milvus async client.
        collection_name: Name of the Milvus collection to store embeddings.
        graph_db_batch_size: Batch size for Neo4j queries.
        embed_batch_size: Batch size for embedder (smaller than molecules due
            to larger model memory requirements).
        model_name: HuggingFace model identifier for ESM2.
        dtype: Data type for model computation and embeddings.
        max_length: Maximum token length for ESM2. Sequences exceeding this
            are skipped by the embedder.
        device: Device for embedder (-1 for all GPUs, 0+ for specific GPU).
        show_progress: Whether to show a progress bar.
        ids: List of protein IDs to embed. If None, all embeddable proteins
            will be embedded.
        max_seq_length: Optional maximum amino acid sequence length filter
            for Neo4j queries (separate from token max_length).

    Returns:
        Statistics dict with keys:
            - fetched: Number of proteins fetched from Neo4j
            - embedded: Number of proteins successfully embedded
            - inserted: Number of embeddings inserted into Milvus
    """
    stats = {"fetched": 0, "embedded": 0, "inserted": 0}

    # Get IDs to embed: use provided list or fetch all embeddable IDs
    if ids is None:
        ids = await fetch_all_embeddable_protein_ids(driver, max_seq_length)
        logger.info(f"Fetched {len(ids)} embeddable protein IDs from Neo4j")
    else:
        logger.info(f"Using {len(ids)} provided protein IDs")

    if not ids:
        logger.warning("No protein IDs to embed")
        return stats

    # Check if collection exists
    existing_collections = await client.list_collections()
    collection_exists = collection_name in existing_collections
    collection_initialized = collection_exists

    # Get total count for progress bar
    total_proteins = len(ids) if show_progress else None

    async with ESM2Embedder(
        model_name=model_name,
        dtype=dtype,
        max_length=max_length,
        device=device,
    ) as embedder:
        # Build pipeline with ID filtering
        sequence_batches = fetch_sequence_batches(
            driver, batch_size=graph_db_batch_size, ids=ids, max_seq_length=max_seq_length
        )
        fetched_counter = [0]
        embedder_batches = rebatch_for_embedder(
            sequence_batches, batch_size=embed_batch_size, fetched_counter=fetched_counter
        )

        # Create progress context
        progress = create_progress() if show_progress else None

        async def process_batches() -> None:
            nonlocal collection_initialized

            task_id = None
            if progress:
                task_id = progress.add_task("Embedding proteins", total=total_proteins)

            async for batch_results in embedder.embed_stream(embedder_batches):
                if not batch_results:
                    continue

                batch_count = len(batch_results)
                stats["embedded"] += batch_count

                # Convert to Milvus records
                records = [
                    {"id": protein_id, "embedding": embedding}
                    for protein_id, embedding in batch_results
                ]

                # Initialize collection on first batch if needed
                if not collection_initialized:
                    await initialize_collection_from_dict(
                        collection_name=collection_name,
                        primary_field_name="id",
                        record=records[0],
                        client=client,
                    )
                    collection_initialized = True
                    logger.info(f"Created collection '{collection_name}'")

                # Insert records
                await insert(
                    collection_name=collection_name,
                    records=records,
                    client=client,
                )
                stats["inserted"] += len(records)

                # Update progress
                if progress and task_id is not None:
                    progress.update(task_id, advance=batch_count)

                logger.debug(
                    f"Inserted {len(records)} protein embeddings",
                    extra={"total_inserted": stats["inserted"]},
                )

        if progress:
            with progress:
                await process_batches()
        else:
            await process_batches()

    # Update fetched count from counter
    stats["fetched"] = fetched_counter[0]

    logger.info(
        f"Embedding complete: {stats['fetched']} fetched, "
        f"{stats['embedded']} embedded, {stats['inserted']} inserted",
        extra=stats,
    )

    return stats


if __name__ == "__main__":
    import asyncio

    from rich import print

    from pyeed.db.milvus import get_async_milvus_client
    from pyeed.db.neo4j import get_async_driver

    async def main() -> None:
        client = get_async_milvus_client()

        async with get_async_driver() as driver:
            # Get proteins not yet embedded
            embeddable_protein_ids = await get_embeddable_proteins(
                driver,
                client,
                "proteins",
                "id",
            )
            print(f"Found {len(embeddable_protein_ids)} proteins to embed")

            if not embeddable_protein_ids:
                print("No proteins to embed")
                await client.close()
                return

            # Embed proteins
            stats = await embed_proteins(
                ids=embeddable_protein_ids,
                driver=driver,
                client=client,
                collection_name="proteins",
                graph_db_batch_size=1000,
                embed_batch_size=32,
                dtype="float16",
                max_length=1024,
                show_progress=True,
                device=0,
            )
            print(f"Embedding complete: {stats}")

        await client.close()

    asyncio.run(main())
