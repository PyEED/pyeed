"""Reaction enrichment stage for protein data."""

from __future__ import annotations

import asyncio
from collections import defaultdict
from typing import Any

from loguru import logger
from rich.progress import Progress, TaskID

from pyeed.db.neo4j import GraphDB
from pyeed.db.queries import (
    _upsert_nodes_with_session,
    create_relationships_batch,
    query_existing_nodes_by_ids,
    remove_list_property_values,
)
from pyeed.ingest.core.pipeline import PipelineRecord
from pyeed.ingest.core.protocol import SENTINEL, PipelineContext
from pyeed.ingest.model import Reaction
from pyeed.ingest.model.protein import Protein
from pyeed.ingest.sources.rhea import RheaClient


class ReactionEnrichmentStage:
    """Batches proteins and enriches them with reaction data.

    Accumulates proteins, extracts unique reaction_ids, fetches reaction data
    from Rhea API, creates Reaction nodes with CATALYZES relationships, and
    forwards Reaction objects to MoleculeStage for molecule enrichment.
    """

    def __init__(
        self,
        db: GraphDB,
        db_semaphore: asyncio.Semaphore,
        batch_size: int = 1000,
        max_concurrent: int = 50,
    ):
        """Initialize reaction enrichment stage.

        Args:
            db: GraphDB instance for Neo4j operations
            db_semaphore: Shared semaphore for DB operations (REQUIRED)
            batch_size: Number of records to accumulate before enriching
            max_concurrent: Maximum concurrent API requests
        """
        self.db = db
        self.batch_size = batch_size
        self.max_concurrent = max_concurrent
        self.rhea_client = RheaClient()
        self._db_semaphore = db_semaphore
        self._progress_lock = asyncio.Lock()
        # Track unique reaction IDs we have already counted in progress.total
        self._seen_reaction_ids: set[str] = set()
        self._total_unique_reactions: int = 0

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Run reaction enrichment: accumulate → spawn background enrichment tasks.

        Args:
            input_queues: Input queues with PipelineRecords
            output_queues: Output queue for Reaction objects (to MoleculeStage)
            context: Pipeline context for shared state
            progress: Progress instance for tracking
            task_id: Task ID for progress updates
        """
        logger.info("Starting reaction enrichment stage")

        input_queue = next(iter(input_queues.values()))
        output_queue = next(iter(output_queues.values())) if output_queues else None
        batch: list[PipelineRecord[Protein]] = []
        enrichment_tasks: list[asyncio.Task[Any]] = []

        while True:
            item = await input_queue.get()

            if item is SENTINEL:
                # Process remaining batch
                if batch:
                    logger.info(f"Processing final reaction batch of {len(batch)} records")
                    task = asyncio.create_task(
                        self._enrich_batch(batch.copy(), output_queue, progress, task_id)
                    )
                    enrichment_tasks.append(task)
                    batch.clear()

                # Wait for all background enrichment tasks
                if enrichment_tasks:
                    logger.info(f"Waiting for {len(enrichment_tasks)} reaction enrichment tasks...")
                    await asyncio.gather(*enrichment_tasks, return_exceptions=True)
                    logger.info("All reaction enrichment tasks completed")

                # Forward sentinel to output queue
                if output_queue is not None:
                    await output_queue.put(SENTINEL)

                break

            # Accumulate records
            batch.append(item)

            # When batch is full, spawn background enrichment task
            if len(batch) >= self.batch_size:
                logger.info(f"Reaction batch full ({len(batch)} records), spawning task")
                task = asyncio.create_task(
                    self._enrich_batch(batch.copy(), output_queue, progress, task_id)
                )
                enrichment_tasks.append(task)
                batch.clear()

        logger.info("Reaction enrichment stage complete")

    async def _enrich_batch(
        self,
        batch: list[PipelineRecord[Protein]],
        output_queue: asyncio.Queue[Any] | None,
        progress: Progress | None = None,
        task_id: TaskID | None = None,
    ) -> None:
        """Enrich a batch of proteins with reaction data.

        Args:
            batch: List of PipelineRecord objects containing Protein data
            output_queue: Queue to forward Reaction objects to MoleculeStage
            progress: Optional Progress instance for updating progress
            task_id: Optional TaskID for tracking progress
        """
        protein_ids = [record.data.id for record in batch]
        logger.debug(f"Enriching reactions for batch of {len(protein_ids)} proteins")

        try:
            reactions, new_reactions = await self._enrich_reactions_batch(
                protein_ids, progress, task_id
            )
            logger.debug(
                f"Reaction enrichment complete for {len(protein_ids)} proteins, "
                f"forwarding {len(reactions)} reactions"
            )

            # Forward Reaction objects to MoleculeStage
            if output_queue is not None:
                for reaction in reactions:
                    await output_queue.put(reaction)

            # Advance progress by number of new unique reaction IDs processed
            if progress is not None and task_id is not None and new_reactions > 0:
                async with self._progress_lock:
                    progress.advance(task_id, new_reactions)
        except Exception:
            logger.exception("Reaction batch enrichment failed")

    async def _enrich_reactions_batch(
        self,
        protein_ids: list[str],
        progress: Progress | None = None,
        task_id: TaskID | None = None,
    ) -> tuple[list[Reaction], int]:
        """Enrich specific proteins with reaction data and return Reaction objects.

        Args:
            protein_ids: List of protein IDs to enrich
            progress: Optional Progress instance for updating progress
            task_id: Optional TaskID for tracking progress

        Returns:
            Tuple of (list of Reaction objects that were processed, count of new unique reactions)
        """
        # Initialize new_reaction_ids to ensure it's always defined
        new_reaction_ids: list[str] = []

        # Step 1: Query proteins for reaction_ids
        async with self.db.async_driver.session() as session:
            query = """
            MATCH (p:Protein)
            WHERE p.id IN $protein_ids AND p.reaction_ids IS NOT NULL AND size(p.reaction_ids) > 0
            RETURN p.id AS id, p.reaction_ids AS reaction_ids
            """
            result = await session.run(query, protein_ids=protein_ids)

            protein_reaction_map: dict[str, list[str]] = {}
            unique_reaction_ids: list[str] = []

            async for record in result:
                protein_id = record["id"]
                reaction_ids = record["reaction_ids"]
                for reaction_id in reaction_ids:
                    if reaction_id not in protein_reaction_map:
                        protein_reaction_map[reaction_id] = []
                        unique_reaction_ids.append(reaction_id)
                    protein_reaction_map[reaction_id].append(protein_id)

            if not unique_reaction_ids:
                return ([], 0)

            # Count only previously unseen reaction IDs for the progress total
            new_reaction_ids = [
                rid for rid in unique_reaction_ids if rid not in self._seen_reaction_ids
            ]
            if new_reaction_ids:
                self._seen_reaction_ids.update(new_reaction_ids)
                if progress is not None and task_id is not None:
                    async with self._progress_lock:
                        self._total_unique_reactions += len(new_reaction_ids)
                        progress.update(task_id, total=self._total_unique_reactions)

            # Check existing reactions
            existing_reaction_ids = await query_existing_nodes_by_ids(
                session, "Reaction", "id", unique_reaction_ids
            )
            reaction_ids_to_fetch = [
                rid for rid in unique_reaction_ids if rid not in existing_reaction_ids
            ]

        logger.debug(f"Reaction IDs to fetch: {len(reaction_ids_to_fetch)}")

        all_reactions: list[Reaction] = []

        # Step 2: Process existing reactions and get their Reaction objects
        if existing_reaction_ids:
            existing_reactions = await self._link_existing_reactions_batch(
                list(existing_reaction_ids), protein_reaction_map
            )
            all_reactions.extend(existing_reactions)

        # Step 3: Fetch and process new reactions
        if reaction_ids_to_fetch:
            reaction_data: list[tuple[dict[str, Any], dict[str, Any]]] = []
            async for table_row, meta_json in self.rhea_client.fetch_reactions(
                reaction_ids_to_fetch, max_concurrent=self.max_concurrent
            ):
                reaction_data.append((table_row, meta_json))

            if reaction_data:
                new_reactions = await self._process_reactions_batch(
                    reaction_data, protein_reaction_map
                )
                all_reactions.extend(new_reactions)

        return (all_reactions, len(new_reaction_ids))

    async def _link_existing_reactions_batch(
        self, reaction_ids: list[str], protein_reaction_map: dict[str, list[str]]
    ) -> list[Reaction]:
        """Link existing reactions to proteins and return Reaction objects.

        Args:
            reaction_ids: List of existing reaction IDs to link
            protein_reaction_map: Mapping of reaction_id → [protein_ids]

        Returns:
            List of Reaction objects
        """
        if not reaction_ids:
            return []

        # Collect all relationships and removals
        all_source_proteins: list[str] = []
        all_target_reactions: list[str] = []
        removal_map: dict[str, list[str]] = defaultdict(list)

        for reaction_id in reaction_ids:
            protein_ids = protein_reaction_map.get(reaction_id, [])
            if not protein_ids:
                continue
            all_source_proteins.extend(protein_ids)
            all_target_reactions.extend([reaction_id] * len(protein_ids))
            for pid in protein_ids:
                removal_map[pid].append(reaction_id)

        if not all_source_proteins:
            return []

        # Query existing Reaction nodes to get their data
        reaction_objects: list[Reaction] = []
        async with self.db.async_driver.session() as session:
            query = """
            MATCH (r:Reaction)
            WHERE r.id IN $reaction_ids
            RETURN r.id AS id, r.description AS description,
                   r.substrate_ids AS substrate_ids, r.product_ids AS product_ids,
                   r.reversible AS reversible
            """
            result = await session.run(query, reaction_ids=reaction_ids)

            async for record in result:
                reaction = Reaction(
                    id=record["id"],
                    description=record.get("description"),
                    substrate_ids=record.get("substrate_ids") or [],
                    product_ids=record.get("product_ids") or [],
                    reversible=record.get("reversible", False),
                )
                reaction_objects.append(reaction)

        # Transaction 1: Create relationships
        async with self._db_semaphore, self.db.async_driver.session() as session:
            await create_relationships_batch(
                session=session,
                source_label="Protein",
                source_field="id",
                source_values=all_source_proteins,
                target_label="Reaction",
                target_field="id",
                target_values=all_target_reactions,
                relationship_type="CATALYZES",
                direction_to_source=True,
                tx_size=5000,
            )

        # Transaction 2: Remove properties
        if removal_map:
            async with self._db_semaphore, self.db.async_driver.session() as session:
                await remove_list_property_values(
                    session=session,
                    label="Protein",
                    unique_field="id",
                    unique_values=list(removal_map.keys()),
                    list_property="reaction_ids",
                    values_to_remove=removal_map,
                )

        return reaction_objects

    async def _process_reactions_batch(
        self,
        reaction_data: list[tuple[dict[str, Any], dict[str, Any]]],
        protein_reaction_map: dict[str, list[str]],
    ) -> list[Reaction]:
        """Process multiple reaction responses in one transaction.

        Args:
            reaction_data: List of (table_row, meta_json) tuples from Rhea API
            protein_reaction_map: Mapping of reaction_id → [protein_ids]

        Returns:
            List of Reaction objects that were processed
        """
        if not reaction_data:
            return []

        # Extract data outside semaphore
        all_reaction_nodes: list[Reaction] = []
        all_source_proteins: list[str] = []
        all_target_reactions: list[str] = []
        removal_map: dict[str, list[str]] = defaultdict(list)

        for table_row, meta_json in reaction_data:
            try:
                reactions = self.rhea_client._extract_reaction(table_row, meta_json, rhea_id=None)
                if not reactions:
                    continue

                reaction = reactions[0]  # Main reaction
                reaction_id = reaction.id

                protein_ids = protein_reaction_map.get(reaction_id, [])
                if not protein_ids:
                    continue

                all_reaction_nodes.append(reaction)
                all_source_proteins.extend(protein_ids)
                all_target_reactions.extend([reaction.id] * len(protein_ids))
                for pid in protein_ids:
                    removal_map[pid].append(reaction.id)
            except Exception as e:
                logger.warning(f"Failed to extract reaction data: {e}", exc_info=True)

        if not all_reaction_nodes:
            return []

        # Single transaction for all reactions
        try:
            async with self._db_semaphore, self.db.async_driver.session() as session:
                # Upsert all reaction nodes
                await _upsert_nodes_with_session(session, all_reaction_nodes)

                # Link all Reactions to Proteins
                if all_source_proteins:
                    await create_relationships_batch(
                        session=session,
                        source_label="Protein",
                        source_field="id",
                        source_values=all_source_proteins,
                        target_label="Reaction",
                        target_field="id",
                        target_values=all_target_reactions,
                        relationship_type="CATALYZES",
                        direction_to_source=True,
                    )

                # Remove all reaction_ids from proteins
                if removal_map:
                    await remove_list_property_values(
                        session=session,
                        label="Protein",
                        unique_field="id",
                        unique_values=list(removal_map.keys()),
                        list_property="reaction_ids",
                        values_to_remove=removal_map,
                    )
        except Exception as e:
            logger.warning(f"Failed to persist reactions batch: {e}", exc_info=True)

        return all_reaction_nodes
