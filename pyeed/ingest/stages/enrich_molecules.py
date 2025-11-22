"""Molecule enrichment stage for reaction data."""

from __future__ import annotations

import asyncio
from typing import Any

from loguru import logger
from rich.progress import Progress, TaskID

from pyeed.db.neo4j import GraphDB
from pyeed.db.queries import (
    _upsert_nodes_with_session,
    create_reaction_molecule_relationships,
    query_existing_nodes_by_ids,
)
from pyeed.ingest.core.protocol import SENTINEL, PipelineContext
from pyeed.ingest.model import Molecule, Reaction
from pyeed.ingest.sources.chebi import ChebiClient


class MoleculeEnrichmentStage:
    """Batches reactions and enriches them with molecule data.

    Receives Reaction objects from ReactionStage, extracts unique molecule IDs
    from substrate_ids and product_ids, fetches molecule data from ChEBI API,
    and creates Molecule nodes with HAS_SUBSTRATE/HAS_PRODUCT relationships.
    """

    def __init__(
        self,
        db: GraphDB,
        db_semaphore: asyncio.Semaphore,
        batch_size: int = 1000,
        max_concurrent: int = 50,
    ):
        """Initialize molecule enrichment stage.

        Args:
            db: GraphDB instance for Neo4j operations
            db_semaphore: Shared semaphore for DB operations (REQUIRED)
            batch_size: Number of reactions to accumulate before enriching
            max_concurrent: Maximum concurrent API requests
        """
        self.db = db
        self.batch_size = batch_size
        self.max_concurrent = max_concurrent
        self.chebi_client = ChebiClient()
        self._db_semaphore = db_semaphore
        self._progress_lock = asyncio.Lock()

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Run molecule enrichment: accumulate → spawn background enrichment tasks.

        Args:
            input_queues: Input queues with Reaction objects
            output_queues: Output queues (not used - terminal stage)
            context: Pipeline context for shared state
            progress: Progress instance for tracking
            task_id: Task ID for progress updates
        """
        logger.info("Starting molecule enrichment stage")

        input_queue = next(iter(input_queues.values()))
        batch: list[Reaction] = []
        enrichment_tasks: list[asyncio.Task[Any]] = []
        cumulative_reaction_count = 0

        while True:
            item = await input_queue.get()

            if item is SENTINEL:
                # Process remaining batch
                if batch:
                    logger.info(f"Processing final molecule batch of {len(batch)} reactions")
                    task = asyncio.create_task(self._enrich_batch(batch.copy(), progress, task_id))
                    enrichment_tasks.append(task)
                    batch.clear()

                # Wait for all background enrichment tasks
                if enrichment_tasks:
                    logger.info(f"Waiting for {len(enrichment_tasks)} molecule enrichment tasks...")
                    await asyncio.gather(*enrichment_tasks, return_exceptions=True)
                    logger.info("All molecule enrichment tasks completed")

                break

            # Accumulate reactions
            if isinstance(item, Reaction):
                batch.append(item)

                # Update cumulative count and progress total dynamically
                if progress is not None and task_id is not None:
                    cumulative_reaction_count += 1
                    async with self._progress_lock:
                        progress.update(task_id, total=cumulative_reaction_count)

                # When batch is full, spawn background enrichment task
                if len(batch) >= self.batch_size:
                    logger.info(f"Molecule batch full ({len(batch)} reactions), spawning task")
                    task = asyncio.create_task(self._enrich_batch(batch.copy(), progress, task_id))
                    enrichment_tasks.append(task)
                    batch.clear()

        logger.info("Molecule enrichment stage complete")

    async def _enrich_batch(
        self,
        batch: list[Reaction],
        progress: Progress | None = None,
        task_id: TaskID | None = None,
    ) -> None:
        """Enrich a batch of reactions with molecule data.

        Args:
            batch: List of Reaction objects
            progress: Optional Progress instance for updating progress
            task_id: Optional TaskID for tracking progress
        """
        logger.debug(f"Enriching molecules for batch of {len(batch)} reactions")

        try:
            await self._enrich_molecules_batch(batch)
            logger.debug(f"Molecule enrichment complete for {len(batch)} reactions")

            # Advance progress by number of reactions processed
            if progress is not None and task_id is not None:
                async with self._progress_lock:
                    progress.advance(task_id, len(batch))
        except Exception:
            logger.exception("Molecule batch enrichment failed")

    async def _enrich_molecules_batch(self, reactions: list[Reaction]) -> None:
        """Enrich reactions with molecule data.

        Args:
            reactions: List of Reaction objects to enrich with molecules
        """
        # Extract all unique molecule IDs from reactions
        substrate_map: dict[str, list[str]] = {}
        product_map: dict[str, list[str]] = {}
        all_chebi_ids: set[str] = set()

        for reaction in reactions:
            if reaction.substrate_ids:
                substrate_map[reaction.id] = reaction.substrate_ids
                all_chebi_ids.update(reaction.substrate_ids)
            if reaction.product_ids:
                product_map[reaction.id] = reaction.product_ids
                all_chebi_ids.update(reaction.product_ids)

        if not all_chebi_ids:
            return

        logger.debug(f"Enriching {len(all_chebi_ids)} unique molecules")

        # Check which molecules already exist
        async with self._db_semaphore, self.db.async_driver.session() as session:
            existing_molecule_ids = await query_existing_nodes_by_ids(
                session, "Molecule", "id", list(all_chebi_ids)
            )

        molecule_ids_to_fetch = [mid for mid in all_chebi_ids if mid not in existing_molecule_ids]

        # Fetch new molecules
        all_molecule_nodes: list[Molecule] = []
        if molecule_ids_to_fetch:
            logger.debug(f"Fetching {len(molecule_ids_to_fetch)} new molecules from ChEBI")
            async for chebi_entry in self.chebi_client.fetch_molecules(
                molecule_ids_to_fetch, batch_size=50
            ):
                try:
                    molecules = self.chebi_client._extract_molecules(chebi_entry)
                    all_molecule_nodes.extend(molecules)
                except Exception as e:
                    logger.warning(f"Failed to extract molecule: {e}", exc_info=True)

        # Split into 2 separate transactions to reduce lock contention
        try:
            # Transaction 1: Upsert molecules
            if all_molecule_nodes:
                async with self._db_semaphore, self.db.async_driver.session() as session:
                    await _upsert_nodes_with_session(session, all_molecule_nodes)

            # Transaction 2: Create all reaction-molecule relationships
            if substrate_map or product_map:
                async with self._db_semaphore, self.db.async_driver.session() as session:
                    await create_reaction_molecule_relationships(
                        session=session,
                        substrate_map=substrate_map,
                        product_map=product_map,
                    )
        except Exception:
            logger.exception("Failed to persist molecules batch")
