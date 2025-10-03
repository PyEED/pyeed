# pyeed/main_enrichment.py
"""
New enrichment method to add to main.py.
This enriches existing Reaction nodes with Molecule relationships.
"""

from __future__ import annotations

import asyncio
import logging

from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TaskID,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from .chebi import ChebiClient
from .database import Database
from .model import MODEL_CLASSES, Molecule
from .rhea import RheaClient

logger = logging.getLogger(__name__)


async def enrich_reactions_with_molecules(
    db: Database,
    concurrency: int = 8,
    save_batch: int = 200,
    sync_schema: bool = False,
) -> None:
    """
    Enrich existing Reaction nodes with Molecule edges via RheaClient and ChebiClient.

    This method:
    1. Queries database for reactions without molecule relationships
    2. Fetches reaction data from Rhea API to extract ChEBI IDs
    3. Batches unique ChEBI IDs and fetches molecule data from ChEBI API
    4. Creates/matches Molecule nodes in the database
    5. Creates HAS_SUBSTRATE and HAS_PRODUCT edges between reactions and molecules

    Args:
        db: Database instance
        concurrency: Max concurrent HTTP requests
        save_batch: Batch size for database saves
        sync_schema: Whether to sync database schema first
    """
    if sync_schema:
        await db.sync_schema(MODEL_CLASSES)

    rhea_client = RheaClient(user_agent="pyeed/1.0")
    chebi_client = ChebiClient(user_agent="pyeed/1.0")

    # Find reactions with RHEA IDs but no molecule relationships
    pquery = """
    MATCH (r:Reaction)
    WHERE r.rhea_id IS NOT NULL
      AND NOT (r)--(:Molecule)
    RETURN r.rhea_id AS rhea_id
    """

    # Count candidates for progress tracking
    count_result = db.query(
        """
        MATCH (r:Reaction)
        WHERE r.rhea_id IS NOT NULL
          AND NOT (r)--(:Molecule)
        RETURN count(r) AS n
    """
    )
    total_candidates = count_result[0]["n"] if count_result else 0
    logger.info(f"Found {total_candidates} reactions to enrich with molecules")

    if total_candidates == 0:
        logger.info("No reactions need enrichment.")
        return

    # Data structures for producer-consumer
    # Queue holds (rhea_id, substrate_chebi_ids, product_chebi_ids)
    queue: asyncio.Queue[tuple[str, list[str], list[str]] | None] = asyncio.Queue(
        maxsize=concurrency * 2
    )
    sem = asyncio.Semaphore(concurrency)

    async def fetch_chebi_ids_for_reaction(
        rhea_id: str,
    ) -> tuple[str, list[str], list[str]] | None:
        """Fetch ChEBI IDs for a reaction from Rhea API."""
        async with sem:
            try:
                # Use rhea.py to get the reaction data
                reaction = await rhea_client.get_reaction(rhea_id)
                if not reaction:
                    logger.warning(f"No reaction data found for RHEA ID: {rhea_id}")
                    return None

                # Extract ChEBI IDs from substrates and products
                substrate_chebi_ids = [mol.chebi_id for mol in reaction.substrates if mol.chebi_id]
                product_chebi_ids = [mol.chebi_id for mol in reaction.products if mol.chebi_id]

                if not substrate_chebi_ids and not product_chebi_ids:
                    logger.warning(f"No ChEBI IDs found for reaction {rhea_id}")
                    return None

                logger.debug(
                    f"Fetched {len(substrate_chebi_ids)} substrates and {len(product_chebi_ids)} products for {rhea_id}"
                )

                return rhea_id, substrate_chebi_ids, product_chebi_ids

            except Exception as e:
                logger.error(f"Failed to fetch reaction {rhea_id}: {e}")
                return None

    async def producer(fetch_task_id: TaskID) -> None:
        """Producer: fetch ChEBI IDs from Rhea API for each reaction."""
        tasks = []
        async for rhea_id in db.async_value_iter(pquery, "rhea_id"):
            tasks.append(asyncio.create_task(fetch_chebi_ids_for_reaction(rhea_id)))

            # Drain in waves to cap memory
            if len(tasks) >= concurrency * 8:
                for coro in asyncio.as_completed(tasks):
                    result = await coro
                    if result:
                        await queue.put(result)
                        progress.update(fetch_task_id, advance=1)
                tasks.clear()

        # Flush remaining tasks
        for coro in asyncio.as_completed(tasks):
            result = await coro
            if result:
                await queue.put(result)
                progress.update(fetch_task_id, advance=1)

        await queue.put(None)  # Signal completion

    async def consumer(save_task_id: TaskID) -> None:
        """Consumer: batch fetch molecules from ChEBI and create edges."""
        # Collect data for batch processing
        reaction_molecule_map: dict[str, tuple[list[str], list[str]]] = {}
        all_chebi_ids: set[str] = set()

        # Phase 1: Collect all data from queue
        while True:
            item = await queue.get()
            if item is None:
                break

            rhea_id, substrate_chebi_ids, product_chebi_ids = item
            reaction_molecule_map[rhea_id] = (substrate_chebi_ids, product_chebi_ids)
            all_chebi_ids.update(substrate_chebi_ids)
            all_chebi_ids.update(product_chebi_ids)

        logger.info(
            f"Collected {len(all_chebi_ids)} unique ChEBI IDs from {len(reaction_molecule_map)} reactions"
        )

        # Phase 2: Batch fetch molecules and save
        chebi_to_molecule = await _fetch_and_save_molecules(
            chebi_client, db, list(all_chebi_ids), save_batch
        )

        # Phase 3: Create edges between reactions and molecules
        await _create_reaction_molecule_edges(
            db, reaction_molecule_map, chebi_to_molecule, save_batch
        )

        # Update progress
        progress.update(save_task_id, advance=len(reaction_molecule_map))

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        transient=False,
    )

    with progress:
        fetch_task = progress.add_task("Fetch Rhea Data", total=total_candidates)
        save_task = progress.add_task("Save Molecules & Edges", total=total_candidates)

        await asyncio.gather(
            producer(fetch_task),
            consumer(save_task),
        )


async def _fetch_and_save_molecules(
    chebi_client: ChebiClient,
    db: Database,
    unique_chebi_list: list[str],
    save_batch: int,
) -> dict[str, Molecule]:
    """Fetch molecules from ChEBI and save to database."""
    chebi_to_molecule: dict[str, Molecule] = {}

    # Fetch molecules in batches
    for i in range(0, len(unique_chebi_list), save_batch):
        batch_ids = unique_chebi_list[i : i + save_batch]
        batch_num = i // save_batch + 1
        total_batches = (len(unique_chebi_list) + save_batch - 1) // save_batch

        try:
            molecules = await chebi_client.get_molecules(batch_ids)
            for mol in molecules:
                chebi_to_molecule[mol.chebi_id] = mol
            logger.info(f"Fetched {len(molecules)} molecules (batch {batch_num}/{total_batches})")
        except Exception as e:
            logger.error(f"Failed to fetch ChEBI batch {i}-{i + save_batch}: {e}")
            # Try individual fetches as fallback
            for chebi_id in batch_ids:
                try:
                    mol = await chebi_client.get_molecule(chebi_id)
                    chebi_to_molecule[chebi_id] = mol
                except Exception as e2:
                    logger.error(f"Failed to fetch individual ChEBI {chebi_id}: {e2}")

    # Save all molecules to database
    if chebi_to_molecule:
        molecules_to_save = list(chebi_to_molecule.values())
        for i in range(0, len(molecules_to_save), save_batch):
            batch = molecules_to_save[i : i + save_batch]
            await db.save_many(batch)
        logger.info(f"Saved {len(molecules_to_save)} molecules to database")

    return chebi_to_molecule


async def _create_reaction_molecule_edges(
    db: Database,
    reaction_molecule_map: dict[str, tuple[list[str], list[str]]],
    chebi_to_molecule: dict[str, Molecule],
    save_batch: int,
) -> None:
    """Create HAS_SUBSTRATE and HAS_PRODUCT edges in batches."""
    substrate_edge_data: list[dict[str, str]] = []
    product_edge_data: list[dict[str, str]] = []

    for rhea_id, (substrate_ids, product_ids) in reaction_molecule_map.items():
        for chebi_id in substrate_ids:
            if chebi_id in chebi_to_molecule:
                substrate_edge_data.append({"rhea_id": rhea_id, "chebi_id": chebi_id})

        for chebi_id in product_ids:
            if chebi_id in chebi_to_molecule:
                product_edge_data.append({"rhea_id": rhea_id, "chebi_id": chebi_id})

    # Batch create edges
    async with db.async_driver.session() as session:
        # Create substrate edges
        for i in range(0, len(substrate_edge_data), save_batch):
            batch = substrate_edge_data[i : i + save_batch]
            await session.run(
                """
                UNWIND $edges AS edge
                MATCH (r:Reaction {rhea_id: edge.rhea_id})
                MATCH (m:Molecule {chebi_id: edge.chebi_id})
                MERGE (r)-[:HAS_SUBSTRATE]->(m)
                """,
                edges=batch,
            )
        logger.info(f"Created {len(substrate_edge_data)} HAS_SUBSTRATE edges")

        # Create product edges
        for i in range(0, len(product_edge_data), save_batch):
            batch = product_edge_data[i : i + save_batch]
            await session.run(
                """
                UNWIND $edges AS edge
                MATCH (r:Reaction {rhea_id: edge.rhea_id})
                MATCH (m:Molecule {chebi_id: edge.chebi_id})
                MERGE (r)-[:HAS_PRODUCT]->(m)
                """,
                edges=batch,
            )
        logger.info(f"Created {len(product_edge_data)} HAS_PRODUCT edges")


# Example usage
if __name__ == "__main__":
    import asyncio

    from pyeed.database import Database

    async def main() -> None:
        db = Database()
        await enrich_reactions_with_molecules(db, sync_schema=True)

    asyncio.run(main())
