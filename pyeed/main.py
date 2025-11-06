from __future__ import annotations

import asyncio
import logging
from collections.abc import Iterable

import httpx
from rich.console import Console
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

from .db.neo4j import Database
from .environment import IN_NOTEBOOK
from .model import MODEL_CLASSES, Protein, Reaction
from .rhea import RheaClient
from .uniprot import UniProtAdapter

logger = logging.getLogger(__name__)

CONSOLE = Console(force_jupyter=IN_NOTEBOOK)


async def ingest_interpro(
    db: Database,
    interpro: str,
    chunk_size: int = 30,
    page_size: int = 30,
    batch_size: int = 200,
    sync_schema: bool = True,
) -> None:
    """Ingests proteins from UniProt by InterPro ID.

    Args:
        db: GraphDatabase
        interpro: InterPro ID (e.g., "IPR002133")
        chunk_size: Number of accessions per batch for fetching
        page_size: Results per page for UniProt API
        batch_size: Number of proteins to batch before upserting to database
        sync_schema: Whether to sync the database schema before ingestion
    """
    adapter = UniProtAdapter()

    logger.info(f"Fetching accessions for InterPro ID: {interpro}")
    async with httpx.AsyncClient() as client:
        accessions = await adapter.fetch_accessions_by_interpro(client, interpro)

    logger.info(f"Found {len(accessions)} accessions for {interpro}")

    # Forward to ingest_uniprot
    await ingest_uniprot(
        db=db,
        accessions=accessions,
        chunk_size=chunk_size,
        page_size=page_size,
        batch_size=batch_size,
        sync_schema=sync_schema,
    )


async def ingest_uniprot(
    db: Database,
    accessions: Iterable[str],
    chunk_size: int = 30,
    page_size: int = 30,
    batch_size: int = 200,
    sync_schema: bool = True,
) -> None:
    """Ingests proteins from UniProt by accession IDs.

    Args:
        db: GraphDatabase
        accessions: Iterable[str]
        chunk_size: int = 30
        page_size: int = 30
        batch_size: int = 200
        sync_schema: bool = False
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
        console=CONSOLE,
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


async def enrich_rhea(db: Database, concurrency: int = 8, save_batch: int = 200) -> None:
    """
    Enrich Reaction nodes (missing Molecule edges) via RheaClient.
    Uses bounded concurrency for HTTP and batched MERGE writes.
    """
    rhea_client = RheaClient(user_agent="")  # add contact

    pquery = """
    MATCH (r:Reaction)
    WHERE r.rhea_id IS NOT NULL
      AND NOT (r)--(:Molecule)
    RETURN r.rhea_id AS rhea_id
    """
    # preflight
    n = db.query("""
        MATCH (n) RETURN count(n) AS n;
    """)
    print(f"enrich_rhea: database={2} candidates={n}", flush=True)

    queue: asyncio.Queue[Reaction | None] = asyncio.Queue(maxsize=concurrency * 2)
    sem = asyncio.Semaphore(concurrency)

    async def fetch_one(rid: str) -> Reaction | None:
        async with sem:
            return await rhea_client.get_reaction(rid)

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
