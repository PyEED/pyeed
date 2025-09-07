from __future__ import annotations

import asyncio
import logging
from typing import Iterable, List, Optional

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

from .database import GraphDatabase
from .fetch.uniprotadapter import UniProtAdapter
from .model import MODEL_CLASSES, Protein

logger = logging.getLogger(__name__)


async def ingest_uniprot(
    db: GraphDatabase,
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

    if not sync_schema:
        await db.sync_schema(MODEL_CLASSES)

    # filter out already-present proteins
    query = """
    MATCH (p:Protein)
    WHERE p.accession_id IN $accessions
    RETURN p.accession_id AS accession_id
    """
    async with db.driver.session() as session:
        result = await session.run(query, accessions=accessions)
        existing = set(await result.value("accession_id"))

    to_fetch = [a for a in accessions if a not in existing]
    if not to_fetch:
        return

    queue: asyncio.Queue[Optional[Protein]] = asyncio.Queue(maxsize=batch_size * 2)
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
                client, to_fetch, chunk_size=chunk_size, size_per_page=page_size
            ):
                prot = adapter.map(rec)
                await queue.put(prot)
                progress.update(fetch_task_id, advance=1)
        # signal done
        await queue.put(None)

    async def consumer(upsert_task_id: TaskID) -> None:
        batch: List[Protein] = []
        while True:
            item = await queue.get()
            if item is None:
                if batch:
                    await db.upsert_nodes(batch)
                    progress.update(upsert_task_id, advance=len(batch))
                break
            batch.append(item)
            if len(batch) >= batch_size:
                await db.upsert_nodes(batch)
                progress.update(upsert_task_id, advance=len(batch))
                batch.clear()

    with progress:
        fetch_task = progress.add_task("Fetch/Map", total=len(to_fetch))
        upsert_task = progress.add_task("Upsert", total=len(to_fetch))

        await asyncio.gather(
            producer(fetch_task),
            consumer(upsert_task),
        )


if __name__ == "__main__":
    import asyncio

    # load accessions from ids.tsv (3rd column or whole line fallback)
    ids: List[str] = []
    path = "ids.tsv"
    with open(path, "r") as f:
        next(f, None)
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split("\t")
            ids.append(parts[2] if len(parts) > 2 else parts[0])

    print(f"Ingesting {len(ids)} proteins")

    async def main() -> None:
        db = GraphDatabase()
        await db.verify_connection()
        try:
            await ingest_uniprot(
                db,
                ids,
                chunk_size=30,
                page_size=30,
                batch_size=200,
            )
        finally:
            await db.close()

    asyncio.run(main())
