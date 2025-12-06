from __future__ import annotations

from collections.abc import AsyncIterator

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
from pyeed.embed.molecule import MoleculeEmbedder


def create_progress() -> Progress:
    """Create a Rich progress bar for molecule embedding."""
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


async def fetch_smiles_batches(
    driver: AsyncDriver,
    batch_size: int = 1000,
    ids: list[str] | None = None,
) -> AsyncIterator[list[tuple[str, str]]]:
    """Fetch molecules with SMILES from Neo4j in batches.

    Queries molecules that have a non-null SMILES attribute, filtering out
    molecules with wildcard notation (asterisks in SMILES). Optionally filters
    by a specific list of molecule IDs.

    Args:
        driver: Neo4j async driver.
        batch_size: Number of molecules to fetch per query.
        ids: Optional list of molecule IDs to filter by. If None, fetches all
            embeddable molecules.

    Yields:
        Batches of (id, smiles) tuples for valid molecules.
    """
    if ids:
        # Query with ID filter
        query = """
        MATCH (m:Molecule)
        WHERE m.id IN $ids
          AND m.smiles IS NOT NULL
          AND NOT m.smiles CONTAINS '*'
        RETURN m.id AS id, m.smiles AS smiles
        SKIP $offset LIMIT $limit
        """
    else:
        # Query all embeddable molecules
        query = """
        MATCH (m:Molecule)
        WHERE m.smiles IS NOT NULL
          AND NOT m.smiles CONTAINS '*'
        RETURN m.id AS id, m.smiles AS smiles
        SKIP $offset LIMIT $limit
        """

    offset = 0
    while True:
        async with driver.session() as session:
            params = {"offset": offset, "limit": batch_size}
            if ids:
                params["ids"] = ids
            result = await session.run(query, params)
            records = await result.data()

        if not records:
            break

        batch = [(r["id"], r["smiles"]) for r in records]
        logger.debug(f"Fetched {len(batch)} molecules from Neo4j", extra={"offset": offset})
        yield batch

        if len(records) < batch_size:
            break

        offset += batch_size


async def rebatch_for_embedder(
    source: AsyncIterator[list[tuple[str, str]]],
    batch_size: int = 128,
    fetched_counter: list[int] | None = None,
) -> AsyncIterator[tuple[list[str], list[str]]]:
    """Rebatch molecules into fixed-size batches for the embedder.

    Takes variable-sized batches from source and yields fixed-size
    (ids, smiles) tuples suitable for the embedder.

    Args:
        source: Async iterator yielding batches of (id, smiles) tuples.
        batch_size: Target batch size for embedder.
        fetched_counter: Optional mutable list to track total fetched count.

    Yields:
        Tuples of (ids_list, smiles_list) with up to batch_size items.
    """
    buffer_ids: list[str] = []
    buffer_smiles: list[str] = []

    async for batch in source:
        for mol_id, smiles in batch:
            buffer_ids.append(mol_id)
            buffer_smiles.append(smiles)
            if fetched_counter is not None:
                fetched_counter[0] += 1

            if len(buffer_ids) >= batch_size:
                yield buffer_ids[:batch_size], buffer_smiles[:batch_size]
                buffer_ids = buffer_ids[batch_size:]
                buffer_smiles = buffer_smiles[batch_size:]

    # Yield remainder
    if buffer_ids:
        yield buffer_ids, buffer_smiles


async def embed_molecules(
    driver: AsyncDriver,
    client: AsyncMilvusClient,
    collection_name: str,
    graph_db_batch_size: int = 1000,
    embed_batch_size: int = 128,
    device: int = -1,
    show_progress: bool = True,
    ids: list[str] | None = None,
) -> dict[str, int]:
    """Embed molecules from Neo4j and store embeddings in Milvus.

    Fetches molecules with SMILES from Neo4j, embeds them using the
    MoleculeEmbedder, and inserts the embeddings into Milvus. Creates
    the collection automatically if it doesn't exist.

    Args:
        driver: Neo4j async driver.
        client: Milvus async client.
        collection_name: Name of the Milvus collection to store embeddings.
        graph_db_batch_size: Batch size for Neo4j queries.
        embed_batch_size: Batch size for embedder.
        device: Device for embedder (-1 for all GPUs, 0+ for specific GPU).
        show_progress: Whether to show a progress bar.
        ids: List of molecule IDs to embed. If None, all embeddable molecules will be embedded.
    Returns:
        Statistics dict with keys:
            - fetched: Number of molecules fetched from Neo4j
            - embedded: Number of molecules successfully embedded
            - inserted: Number of embeddings inserted into Milvus
    """
    stats = {"fetched": 0, "embedded": 0, "inserted": 0}

    # Get IDs to embed: use provided list or fetch all embeddable IDs
    if ids is None:
        ids = await fetch_all_embeddable_ids(driver)
        logger.info(f"Fetched {len(ids)} embeddable molecule IDs from Neo4j")
    else:
        logger.info(f"Using {len(ids)} provided molecule IDs")

    if not ids:
        logger.warning("No molecule IDs to embed")
        return stats

    # Check if collection exists
    existing_collections = await client.list_collections()
    collection_exists = collection_name in existing_collections
    collection_initialized = collection_exists

    # Get total count for progress bar
    total_molecules = len(ids) if show_progress else None

    async with MoleculeEmbedder(device=device) as embedder:
        # Build pipeline with ID filtering
        smiles_batches = fetch_smiles_batches(driver, batch_size=graph_db_batch_size, ids=ids)
        fetched_counter = [0]
        embedder_batches = rebatch_for_embedder(
            smiles_batches, batch_size=embed_batch_size, fetched_counter=fetched_counter
        )

        # Create progress context
        progress = create_progress() if show_progress else None

        async def process_batches() -> None:
            nonlocal collection_initialized

            task_id = None
            if progress:
                task_id = progress.add_task("Embedding molecules", total=total_molecules)

            async for batch_results in embedder.embed_stream(embedder_batches):
                if not batch_results:
                    continue

                batch_count = len(batch_results)
                stats["embedded"] += batch_count

                # Convert to Milvus records
                records = [
                    {"id": mol_id, "embedding": embedding} for mol_id, embedding in batch_results
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
                    f"Inserted {len(records)} embeddings",
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


async def fetch_all_embeddable_ids(driver: AsyncDriver) -> list[str]:
    """Fetch all molecule IDs with valid SMILES from Neo4j.

    Args:
        driver: Neo4j async driver.

    Returns:
        List of molecule IDs with non-null SMILES without wildcards.
    """
    query = """
    MATCH (m:Molecule)
    WHERE m.smiles IS NOT NULL
      AND NOT m.smiles CONTAINS '*'
    RETURN m.id AS id
    """

    ids: list[str] = []
    async with driver.session() as session:
        result = await session.run(query)
        async for record in result:
            if record["id"] is not None:
                ids.append(record["id"])

    return ids


async def count_embeddable_molecules(driver: AsyncDriver, ids: list[str] | None = None) -> int:
    """Count molecules with valid SMILES for embedding.

    Args:
        driver: Neo4j async driver.
        ids: Optional list of molecule IDs to count. If None, counts all
            embeddable molecules.

    Returns:
        Total count of molecules with non-null SMILES without wildcards.
    """
    if ids:
        query = """
        MATCH (m:Molecule)
        WHERE m.id IN $ids
          AND m.smiles IS NOT NULL
          AND NOT m.smiles CONTAINS '*'
        RETURN count(m) AS count
        """
    else:
        query = """
        MATCH (m:Molecule)
        WHERE m.smiles IS NOT NULL
          AND NOT m.smiles CONTAINS '*'
        RETURN count(m) AS count
        """

    async with driver.session() as session:
        params = {"ids": ids} if ids else {}
        result = await session.run(query, params)
        record = await result.single()
        return record["count"] if record else 0


async def get_embeddable_molecules(
    driver: AsyncDriver,
    client: AsyncMilvusClient,
    collection_name: str,
    primary_field_name: str,
) -> list[str]:
    """Return molecule IDs that need embedding (not yet in Milvus).

    Fetches all molecules with valid SMILES from Neo4j and filters out
    those that already exist in the Milvus collection.

    Args:
        driver: Neo4j async driver.
        client: Milvus async client.
        collection_name: Name of the Milvus collection to check.
        primary_field_name: Name of the primary key field in Milvus.

    Returns:
        List of molecule IDs that need embedding (not in collection).
        Returns all embeddable molecule IDs if collection doesn't exist.
    """
    query = """
    MATCH (m:Molecule)
    WHERE m.smiles IS NOT NULL
      AND NOT m.smiles CONTAINS '*'
    RETURN m.id AS id
    """

    # Fetch all embeddable molecule IDs from Neo4j
    all_ids: list[str] = []
    async with driver.session() as session:
        result = await session.run(query)
        async for record in result:
            if record["id"] is not None:
                all_ids.append(record["id"])

    if not all_ids:
        return []

    # Check if collection exists
    existing_collections = await client.list_collections()
    if collection_name not in existing_collections:
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
        f"Found {len(all_ids)} embeddable molecules in Neo4j, "
        f"{len(existing_ids)} already in Milvus, "
        f"{len(all_ids) - len(existing_ids)} need embedding",
        extra={
            "neo4j_count": len(all_ids),
            "milvus_count": len(existing_ids),
            "to_embed": len(all_ids) - len(existing_ids),
        },
    )

    # Return IDs that need embedding
    return [mol_id for mol_id in all_ids if mol_id not in existing_ids]


if __name__ == "__main__":
    import asyncio

    from rich import print

    from pyeed.db.milvus import get_async_milvus_client
    from pyeed.db.neo4j import get_async_driver

    async def main() -> None:
        client = get_async_milvus_client()

        async with get_async_driver() as driver:
            # Count embeddable molecules
            embeddable_mol_ids = await get_embeddable_molecules(
                driver,
                client,
                "molecules",
                "id",
            )
            print(f"Found {len(embeddable_mol_ids)} embeddable molecules")

            # Embed molecules
            collection_name = "molecules"

            stats = await embed_molecules(
                ids=embeddable_mol_ids,
                driver=driver,
                client=client,
                collection_name=collection_name,
                graph_db_batch_size=1280,
                embed_batch_size=128,
                show_progress=True,
                device=2,
            )
        await client.close()
        logger.info(f"Molecule embedding complete: {stats}")

    asyncio.run(main())
