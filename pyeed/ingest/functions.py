from __future__ import annotations

import asyncio
import logging
from collections.abc import Iterable

import httpx
from rich.progress import Progress, TaskID

from ..db.neo4j import Database
from ..model import MODEL_CLASSES, Protein, Reaction
from .chebi import ChebiClient
from .progress import create_progress
from .rhea import RheaClient
from .uniprot import UniProtAdapter

__all__ = [
    "enrich_molecules",
    "enrich_reactions",
    "fetch_proteins_by_ids",
    "fetch_proteins_by_interpro",
    "ingest_full_pipeline",
]

logger = logging.getLogger(__name__)


async def fetch_proteins_by_ids(
    db: Database,
    accessions: Iterable[str],
    *,
    chunk_size: int = 30,
    page_size: int = 30,
    batch_size: int = 200,
    show_progress: bool = True,
    skip_schema_sync: bool = False,
) -> list[Protein]:
    """Fetch proteins from UniProt by accession IDs.

    Args:
        db: Database instance to save proteins to.
        accessions: UniProt accession IDs to fetch.
        chunk_size: Number of accessions per API request batch.
        page_size: Results per page for UniProt API.
        batch_size: Number of proteins to batch before database upsert.
        show_progress: Display progress bars.
        skip_schema_sync: Skip automatic schema synchronization.

    Returns:
        List of fetched Protein instances.

    Example:
        >>> db = Database()
        >>> proteins = await fetch_proteins_by_ids(
        ...     db,
        ...     ["P12345", "Q9Y6K9"],
        ... )
        >>> len(proteins)
        2
    """
    accessions = list(accessions)

    if not skip_schema_sync:
        await db.sync_schema(MODEL_CLASSES)

    queue: asyncio.Queue[Protein | None] = asyncio.Queue(maxsize=batch_size * 2)
    adapter = UniProtAdapter()
    proteins: list[Protein] = []

    async def producer(fetch_task_id: TaskID | None) -> None:
        async with httpx.AsyncClient() as client:
            async for rec in adapter.fetch_accessions(
                client, accessions, chunk_size=chunk_size, size_per_page=page_size
            ):
                prot = adapter.map(rec)
                proteins.append(prot)
                await queue.put(prot)
                if fetch_task_id is not None:
                    progress.update(fetch_task_id, advance=1)
        await queue.put(None)

    async def consumer(upsert_task_id: TaskID | None) -> None:
        batch: list[Protein] = []
        while True:
            item = await queue.get()
            if item is None:
                if batch:
                    await db.save_many(batch)
                    if upsert_task_id is not None:
                        progress.update(upsert_task_id, advance=len(batch))
                break
            batch.append(item)
            if len(batch) >= batch_size:
                await db.save_many(batch)
                if upsert_task_id is not None:
                    progress.update(upsert_task_id, advance=len(batch))
                batch.clear()

    if show_progress:
        progress = create_progress()
        with progress:
            fetch_task = progress.add_task("Fetch/Map", total=len(accessions))
            upsert_task = progress.add_task("Upsert", total=len(accessions))
            await asyncio.gather(producer(fetch_task), consumer(upsert_task))
    else:
        await asyncio.gather(producer(None), consumer(None))

    logger.info(f"Fetched and saved {len(proteins)} proteins")
    return proteins


async def fetch_proteins_by_interpro(
    db: Database,
    interpro_id: str,
    *,
    chunk_size: int = 30,
    page_size: int = 30,
    batch_size: int = 200,
    show_progress: bool = True,
    skip_schema_sync: bool = False,
) -> list[Protein]:
    """Fetch proteins from UniProt by InterPro ID.

    First queries UniProt SPARQL endpoint to get all accessions linked to the
    InterPro ID, then fetches full protein data for each accession.

    Args:
        db: Database instance to save proteins to.
        interpro_id: InterPro ID (e.g., "IPR002133").
        chunk_size: Number of accessions per API request batch.
        page_size: Results per page for UniProt API.
        batch_size: Number of proteins to batch before database upsert.
        show_progress: Display progress bars.
        skip_schema_sync: Skip automatic schema synchronization.

    Returns:
        List of fetched Protein instances.

    Example:
        >>> db = Database()
        >>> proteins = await fetch_proteins_by_interpro(
        ...     db,
        ...     "IPR002133",
        ... )
        >>> len(proteins) > 0
        True
    """
    adapter = UniProtAdapter()

    logger.info(f"Fetching accessions for InterPro ID: {interpro_id}")
    async with httpx.AsyncClient() as client:
        accessions = await adapter.fetch_accessions_by_interpro(client, interpro_id)

    logger.info(f"Found {len(accessions)} accessions for {interpro_id}")

    return await fetch_proteins_by_ids(
        db,
        accessions,
        chunk_size=chunk_size,
        page_size=page_size,
        batch_size=batch_size,
        show_progress=show_progress,
        skip_schema_sync=skip_schema_sync,
    )


async def enrich_reactions(
    db: Database,
    *,
    concurrency: int = 8,
    batch_size: int = 200,
    show_progress: bool = True,
    progress: Progress | None = None,
) -> int:
    """Enrich Reaction nodes with Molecule edges via Rhea.

    Queries for reactions missing molecule edges, fetches reaction details
    from Rhea, and updates the database with substrate and product molecules.

    Args:
        db: Database instance.
        concurrency: Number of concurrent HTTP requests to Rhea.
        batch_size: Number of reactions to batch before database upsert.
        show_progress: Display progress bars.
        progress: Optional existing Progress instance to share progress context.

    Returns:
        Number of reactions enriched.

    Example:
        >>> db = Database()
        >>> count = await enrich_reactions(db)
        >>> count >= 0
        True
    """
    rhea_client = RheaClient(user_agent="pyeed/1.0")

    # Find reactions missing molecule edges
    pquery = """
    MATCH (r:Reaction)
    WHERE r.rhea_id IS NOT NULL
      AND NOT (r)--(:Molecule)
    RETURN r.rhea_id AS rhea_id
    """

    # Query total count for progress tracking
    count_result = db.query("""
        MATCH (r:Reaction)
        WHERE r.rhea_id IS NOT NULL
          AND NOT (r)--(:Molecule)
        RETURN count(r) AS n
    """)
    total_count = count_result[0]["n"] if count_result else 0

    queue: asyncio.Queue[Reaction | None] = asyncio.Queue(maxsize=concurrency * 2)
    sem = asyncio.Semaphore(concurrency)
    enriched_count = 0

    async def fetch_one(rid: str) -> Reaction | None:
        async with sem:
            return await rhea_client.get_reaction(rid)

    async def producer(fetch_task_id: TaskID | None, prog: Progress | None) -> None:
        tasks = []
        async for rid in db.async_value_iter(pquery, "rhea_id"):
            tasks.append(asyncio.create_task(fetch_one(rid)))
            if len(tasks) >= concurrency * 8:
                for coro in asyncio.as_completed(tasks):
                    rx = await coro
                    if rx:
                        await queue.put(rx)
                        if fetch_task_id is not None and prog is not None:
                            prog.update(fetch_task_id, advance=1)
                tasks.clear()
        for coro in asyncio.as_completed(tasks):
            rx = await coro
            if rx:
                await queue.put(rx)
                if fetch_task_id is not None and prog is not None:
                    prog.update(fetch_task_id, advance=1)
        await queue.put(None)

    async def consumer(save_task_id: TaskID | None, prog: Progress | None) -> None:
        nonlocal enriched_count
        batch: list[Reaction] = []
        while True:
            item = await queue.get()
            if item is None:
                if batch:
                    await db.save_many(batch)
                    enriched_count += len(batch)
                    if save_task_id is not None and prog is not None:
                        prog.update(save_task_id, advance=len(batch))
                break
            batch.append(item)
            if len(batch) >= batch_size:
                await db.save_many(batch)
                enriched_count += len(batch)
                if save_task_id is not None and prog is not None:
                    prog.update(save_task_id, advance=len(batch))
                batch.clear()

    if show_progress:
        progress_instance = create_progress(progress=progress)
        # Only use 'with' context if we created a new progress instance
        if progress is None:
            with progress_instance:
                t_fetch = progress_instance.add_task("Fetch+Map (Rhea)", total=total_count)
                t_save = progress_instance.add_task("Upsert (Neo4j)", total=total_count)
                await asyncio.gather(
                    producer(t_fetch, progress_instance), consumer(t_save, progress_instance)
                )
        else:
            # Shared progress - don't use context manager
            t_fetch = progress_instance.add_task("Fetch+Map (Rhea)", total=total_count)
            t_save = progress_instance.add_task("Upsert (Neo4j)", total=total_count)
            await asyncio.gather(
                producer(t_fetch, progress_instance), consumer(t_save, progress_instance)
            )
    else:
        await asyncio.gather(producer(None, None), consumer(None, None))

    logger.info(f"Enriched {enriched_count} reactions with molecule data")
    return enriched_count


async def enrich_molecules(
    db: Database,
    *,
    batch_size: int = 50,
    show_progress: bool = True,
    progress: Progress | None = None,
) -> int:
    """Enrich Molecule nodes with ChEBI metadata (name, SMILES, InChI).

    Queries for molecules missing metadata, fetches details from ChEBI,
    and updates the database.

    Args:
        db: Database instance.
        batch_size: Number of molecules to fetch per ChEBI API call.
        show_progress: Display progress bars.
        progress: Optional existing Progress instance to share progress context.

    Returns:
        Number of molecules enriched.

    Example:
        >>> db = Database()
        >>> count = await enrich_molecules(db)
        >>> count >= 0
        True
    """
    chebi_client = ChebiClient(user_agent="pyeed/1.0")

    # Find molecules missing metadata (no name, SMILES, or InChI)
    mquery = """
    MATCH (m:Molecule)
    WHERE m.chebi_id IS NOT NULL
      AND (m.name IS NULL OR m.smiles IS NULL OR m.inchi IS NULL)
    RETURN m.chebi_id AS chebi_id
    """

    # Collect all chebi_ids first
    chebi_ids: list[str] = []
    async for cid in db.async_value_iter(mquery, "chebi_id"):
        chebi_ids.append(cid)

    if not chebi_ids:
        logger.info("No molecules need enrichment")
        return 0

    enriched_count = 0

    async def process_batch(
        batch_ids: list[str], task_id: TaskID | None, prog: Progress | None
    ) -> None:
        nonlocal enriched_count
        try:
            molecules = await chebi_client.get_molecules(batch_ids)
            await db.save_many(molecules)
            enriched_count += len(molecules)
            if task_id is not None and prog is not None:
                prog.update(task_id, advance=len(molecules))
        except Exception as e:
            logger.warning(f"Failed to enrich batch of {len(batch_ids)} molecules: {e}")
            # Try individual molecules as fallback
            for chebi_id in batch_ids:
                try:
                    mol = await chebi_client.get_molecule(chebi_id)
                    await db.save_many([mol])
                    enriched_count += 1
                    if task_id is not None and prog is not None:
                        prog.update(task_id, advance=1)
                except Exception as e2:
                    logger.error(f"Failed to enrich molecule {chebi_id}: {e2}")

    if show_progress:
        progress_instance = create_progress(progress=progress)
        # Only use 'with' context if we created a new progress instance
        if progress is None:
            with progress_instance:
                task = progress_instance.add_task("Enrich Molecules (ChEBI)", total=len(chebi_ids))
                batches = [
                    chebi_ids[i : i + batch_size] for i in range(0, len(chebi_ids), batch_size)
                ]
                for batch in batches:
                    await process_batch(batch, task, progress_instance)
        else:
            # Shared progress - don't use context manager
            task = progress_instance.add_task("Enrich Molecules (ChEBI)", total=len(chebi_ids))
            batches = [chebi_ids[i : i + batch_size] for i in range(0, len(chebi_ids), batch_size)]
            for batch in batches:
                await process_batch(batch, task, progress_instance)
    else:
        batches = [chebi_ids[i : i + batch_size] for i in range(0, len(chebi_ids), batch_size)]
        for batch in batches:
            await process_batch(batch, None, None)

    logger.info(f"Enriched {enriched_count} molecules with ChEBI metadata")
    return enriched_count


async def ingest_full_pipeline(
    db: Database,
    accessions: list[str] | None = None,
    interpro_id: str | None = None,
    *,
    include_reactions: bool = True,
    include_molecules: bool = True,
    chunk_size: int = 30,
    page_size: int = 30,
    batch_size: int = 200,
    show_progress: bool = True,
    skip_schema_sync: bool = False,
) -> dict[str, int]:
    """All-in-one ingestion pipeline: fetch proteins + optionally enrich reactions + molecules.

    Provide either accessions or interpro_id (not both).

    Args:
        db: Database instance.
        accessions: List of UniProt accession IDs (mutually exclusive with interpro_id).
        interpro_id: InterPro ID (mutually exclusive with accessions).
        include_reactions: Enrich reactions with Rhea data.
        include_molecules: Enrich molecules with ChEBI data.
        chunk_size: Number of accessions per API request batch.
        page_size: Results per page for UniProt API.
        batch_size: Number of items to batch before database upsert.
        show_progress: Display progress bars.
        skip_schema_sync: Skip automatic schema synchronization.

    Returns:
        Dictionary with counts: {"proteins": int, "reactions": int, "molecules": int}

    Raises:
        ValueError: If neither or both accessions and interpro_id are provided.

    Example:
        >>> db = Database()
        >>> stats = await ingest_full_pipeline(
        ...     db,
        ...     interpro_id="IPR002133",
        ...     include_reactions=True,
        ... )
        >>> stats["proteins"] > 0
        True
    """
    if (accessions is None) == (interpro_id is None):
        raise ValueError("Provide exactly one of: accessions or interpro_id")

    stats: dict[str, int] = {"proteins": 0, "reactions": 0, "molecules": 0}

    # Fetch proteins
    if accessions is not None:
        proteins = await fetch_proteins_by_ids(
            db,
            accessions,
            chunk_size=chunk_size,
            page_size=page_size,
            batch_size=batch_size,
            show_progress=show_progress,
            skip_schema_sync=skip_schema_sync,
        )
    else:
        proteins = await fetch_proteins_by_interpro(
            db,
            interpro_id,  # type: ignore
            chunk_size=chunk_size,
            page_size=page_size,
            batch_size=batch_size,
            show_progress=show_progress,
            skip_schema_sync=skip_schema_sync,
        )

    stats["proteins"] = len(proteins)

    # Create shared progress instance for enrichment steps if needed
    if show_progress and (include_reactions or include_molecules):
        progress = create_progress()
        with progress:
            # Enrich reactions
            if include_reactions:
                stats["reactions"] = await enrich_reactions(
                    db, batch_size=batch_size, show_progress=show_progress, progress=progress
                )

            # Enrich molecules
            if include_molecules:
                stats["molecules"] = await enrich_molecules(
                    db, batch_size=batch_size, show_progress=show_progress, progress=progress
                )
    else:
        # No shared progress needed
        if include_reactions:
            stats["reactions"] = await enrich_reactions(
                db, batch_size=batch_size, show_progress=show_progress
            )

        if include_molecules:
            stats["molecules"] = await enrich_molecules(
                db, batch_size=batch_size, show_progress=show_progress
            )

    logger.info(f"Pipeline complete: {stats}")
    return stats
