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

from .database import Database
from .model import MODEL_CLASSES, Protein
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
