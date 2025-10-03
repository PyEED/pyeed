from __future__ import annotations

import asyncio
import logging
from collections.abc import Iterable

import httpx
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
from .model import MODEL_CLASSES, Protein, Reaction
from .rhea import RheaClient
from .uniprot import UniProtAdapter

logger = logging.getLogger(__name__)


async def ingest_uniprot(
    db: Database,
    accessions: Iterable[str],
    chunk_size: int = 30,
    page_size: int = 30,
    batch_size: int = 200,
    sync_schema: bool = False,
) -> None:
    """Ingests proteins from UniProt by accession IDs.

    Args:
        db: GraphDatabase
        accessions: Iterable[str]
        chunk_size: int = 30
        page_size: int = 30
        batch_size: int = 200
        schema_synced: bool = False
    """
    accessions = list(accessions)

    if sync_schema:
        await db.sync_schema(MODEL_CLASSES)

    queue: asyncio.Queue[Protein | None] = asyncio.Queue(maxsize=batch_size * 2)
    adapter = UniProtAdapter()

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        transient=False,
    )

    async def producer(fetch_task_id: TaskID) -> None:
        async with httpx.AsyncClient() as client:
            async for rec in adapter.fetch_accessions(
                client, accessions, chunk_size=chunk_size, size_per_page=page_size
            ):
                prot = adapter.map(rec)
                await queue.put(prot)
                progress.update(fetch_task_id, advance=1)
        # signal done
        await queue.put(None)

    async def consumer(upsert_task_id: TaskID) -> None:
        batch: list[Protein] = []
        while True:
            item = await queue.get()
            if item is None:
                if batch:
                    await db.save_many(batch)
                    progress.update(upsert_task_id, advance=len(batch))
                break
            batch.append(item)
            if len(batch) >= batch_size:
                await db.save_many(batch)
                progress.update(upsert_task_id, advance=len(batch))
                batch.clear()

    with progress:
        fetch_task = progress.add_task("Fetch/Map", total=len(accessions))
        upsert_task = progress.add_task("Upsert", total=len(accessions))

        await asyncio.gather(
            producer(fetch_task),
            consumer(upsert_task),
        )


async def enrich_rhea(
    db: Database,
    concurrency: int = 8,
    save_batch: int = 200,
) -> None:
    """
    Enrich Reaction nodes (missing Molecule edges) via RheaClient.
    Uses bounded concurrency for HTTP and batched MERGE writes.
    """
    rhea_client = RheaClient(user_agent="pyeed/1.0")
    chebi_client = ChebiClient(user_agent="pyeed/1.0")

    # Cache for ChEBI molecule data to avoid duplicate fetches
    chebi_cache: dict[str, dict[str, str | None] | None] = {}

    pquery = """
    MATCH (r:Reaction)
    WHERE r.rhea_id IS NOT NULL
      AND NOT (r)--(:Molecule)
    RETURN r.rhea_id AS rhea_id
    """

    queue: asyncio.Queue[Reaction | None] = asyncio.Queue(maxsize=concurrency * 2)
    sem = asyncio.Semaphore(concurrency)

    def chebi_enricher(chebi_id: str) -> dict[str, str | None] | None:
        """Synchronous enricher that pulls from cache populated by async fetch."""
        return chebi_cache.get(chebi_id)

    async def fetch_one(rid: str) -> Reaction | None:
        async with sem:
            # First fetch reaction structure without enrichment
            reaction = await rhea_client.get_reaction(rid, chebi_enricher=None)
            if not reaction:
                return None

            # Collect all ChEBI IDs from this reaction
            chebi_ids = [mol.chebi_id for mol in reaction.substrates + reaction.products]

            # Fetch missing ChEBI data
            missing_ids = [cid for cid in chebi_ids if cid not in chebi_cache]
            if missing_ids:
                try:
                    molecules = await chebi_client.get_molecules(missing_ids)
                    for mol in molecules:
                        chebi_cache[mol.chebi_id] = {
                            "name": mol.name,
                            "smiles": mol.smiles,
                            "inchi": mol.inchi,
                        }
                except Exception as e:
                    logger.warning(f"Failed to enrich ChEBI IDs {missing_ids}: {e}")
                    # Cache None for failed IDs to avoid retrying
                    for cid in missing_ids:
                        if cid not in chebi_cache:
                            chebi_cache[cid] = None

            # Now re-fetch reaction with enriched cache
            return await rhea_client.get_reaction(rid, chebi_enricher=chebi_enricher)

    async def producer(fetch_task_id: TaskID) -> None:
        tasks = []
        async for rid in db.async_value_iter(pquery, "rhea_id"):
            tasks.append(asyncio.create_task(fetch_one(rid)))
            # drain in waves to cap memory even if 10k ids
            if len(tasks) >= concurrency * 8:
                for coro in asyncio.as_completed(tasks):
                    rx = await coro
                    if rx:
                        await queue.put(rx)
                        progress.update(fetch_task_id, advance=1)
                tasks.clear()
        # flush remaining
        for coro in asyncio.as_completed(tasks):
            rx = await coro
            if rx:
                await queue.put(rx)
                progress.update(fetch_task_id, advance=1)
        await queue.put(None)

    async def consumer(save_task_id: TaskID) -> None:
        batch: list[Reaction] = []
        while True:
            item = await queue.get()
            if item is None:
                if batch:
                    await db.save_many(batch)  # must MERGE on rhea_id / chebi_id
                    progress.update(save_task_id, advance=len(batch))
                break
            batch.append(item)
            if len(batch) >= save_batch:
                await db.save_many(batch)
                progress.update(save_task_id, advance=len(batch))
                batch.clear()

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        transient=False,
    )

    # unknown total → leave totals None
    with progress:
        t_fetch = progress.add_task("Fetch+Map (Rhea)", total=None)
        t_save = progress.add_task("Upsert (Neo4j)", total=None)
        await asyncio.gather(producer(t_fetch), consumer(t_save))


if __name__ == "__main__":
    import asyncio

    from pyeed.database import Database

    print("Starting")

    db = Database()
    asyncio.run(enrich_rhea(db))
