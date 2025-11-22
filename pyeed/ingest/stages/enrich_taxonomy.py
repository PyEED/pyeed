"""Taxonomy enrichment stage for protein data."""

from __future__ import annotations

import asyncio
from collections import defaultdict
from typing import Any

import httpx
from loguru import logger
from rich.progress import Progress, TaskID

from pyeed.db.neo4j import GraphDB
from pyeed.db.queries import (
    _upsert_nodes_with_session,
    create_relationships_batch,
    create_taxonomy_hierarchy,
    query_existing_nodes_by_ids,
    remove_list_property_values,
)
from pyeed.ingest.core.pipeline import PipelineRecord
from pyeed.ingest.core.protocol import SENTINEL, PipelineContext
from pyeed.ingest.model import Taxon
from pyeed.ingest.model.protein import Protein
from pyeed.ingest.sources.taxonomy import UniProtTaxonomyAdapter


class TaxonomyEnrichmentStage:
    """Batches proteins and enriches them with taxonomy data.

    Accumulates proteins, extracts unique taxon_ids, fetches taxonomy data
    from UniProt, and creates Taxon nodes with ORIGINATES_FROM relationships.
    """

    def __init__(
        self,
        db: GraphDB,
        db_semaphore: asyncio.Semaphore,
        batch_size: int = 1000,
        max_concurrent: int = 50,
    ):
        """Initialize taxonomy enrichment stage.

        Args:
            db: GraphDB instance for Neo4j operations
            db_semaphore: Shared semaphore for DB operations (REQUIRED)
            batch_size: Number of records to accumulate before enriching
            max_concurrent: Maximum concurrent API requests
        """
        self.db = db
        self.batch_size = batch_size
        self.max_concurrent = max_concurrent
        self.taxonomy_adapter = UniProtTaxonomyAdapter()
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
        """Run taxonomy enrichment: accumulate → spawn background enrichment tasks.

        Args:
            input_queues: Input queues with PipelineRecords
            output_queues: Output queues (not used - terminal stage)
            context: Pipeline context for shared state
            progress: Progress instance for tracking
            task_id: Task ID for progress updates
        """
        logger.info("Starting taxonomy enrichment stage")

        input_queue = next(iter(input_queues.values()))
        batch: list[PipelineRecord[Protein]] = []
        enrichment_tasks: list[asyncio.Task[Any]] = []
        cumulative_protein_count = 0

        while True:
            item = await input_queue.get()

            if item is SENTINEL:
                # Process remaining batch
                if batch:
                    logger.info(f"Processing final taxonomy batch of {len(batch)} records")
                    task = asyncio.create_task(self._enrich_batch(batch.copy(), progress, task_id))
                    enrichment_tasks.append(task)
                    batch.clear()

                # Wait for all background enrichment tasks
                if enrichment_tasks:
                    logger.info(f"Waiting for {len(enrichment_tasks)} taxonomy enrichment tasks...")
                    await asyncio.gather(*enrichment_tasks, return_exceptions=True)
                    logger.info("All taxonomy enrichment tasks completed")

                break

            # Accumulate records
            batch.append(item)

            # Update cumulative count and progress total dynamically
            if progress is not None and task_id is not None:
                cumulative_protein_count += 1
                async with self._progress_lock:
                    progress.update(task_id, total=cumulative_protein_count)

            # When batch is full, spawn background enrichment task
            if len(batch) >= self.batch_size:
                logger.info(f"Taxonomy batch full ({len(batch)} records), spawning task")
                task = asyncio.create_task(self._enrich_batch(batch.copy(), progress, task_id))
                enrichment_tasks.append(task)
                batch.clear()

        logger.info("Taxonomy enrichment stage complete")

    async def _enrich_batch(
        self,
        batch: list[PipelineRecord[Protein]],
        progress: Progress | None = None,
        task_id: TaskID | None = None,
    ) -> None:
        """Enrich a batch of proteins with taxonomy data.

        Args:
            batch: List of PipelineRecord objects containing Protein data
            progress: Optional Progress instance for updating progress
            task_id: Optional TaskID for tracking progress
        """
        protein_ids = [record.data.id for record in batch]
        logger.debug(f"Enriching taxonomy for batch of {len(protein_ids)} proteins")

        try:
            await self._enrich_taxonomy_batch(protein_ids)
            logger.debug(f"Taxonomy enrichment complete for {len(protein_ids)} proteins")

            # Advance progress by number of proteins processed
            if progress is not None and task_id is not None:
                async with self._progress_lock:
                    progress.advance(task_id, len(batch))
        except Exception:
            logger.exception("Taxonomy batch enrichment failed")

    async def _enrich_taxonomy_batch(self, protein_ids: list[str]) -> None:
        """Enrich specific proteins with taxonomy data.

        Args:
            protein_ids: List of protein IDs to enrich
        """
        logger.info(f"Starting taxonomy enrichment for {len(protein_ids)} proteins")

        # Step 1: Query proteins for taxon_ids
        async with self.db.async_driver.session() as session:
            query = """
            MATCH (p:Protein)
            WHERE p.id IN $protein_ids AND p.taxon_ids IS NOT NULL AND size(p.taxon_ids) > 0
            RETURN p.id AS id, p.taxon_ids AS taxon_ids
            """
            result = await session.run(query, protein_ids=protein_ids)

            protein_taxon_map: dict[str, list[str]] = defaultdict(list)
            proteins_found = 0
            async for record in result:
                protein_id = record["id"]
                taxon_ids = record["taxon_ids"]
                proteins_found += 1
                for taxon_id in taxon_ids:
                    protein_taxon_map[str(taxon_id)].append(protein_id)

            logger.info(
                f"Found {proteins_found} proteins with taxon_ids out of {len(protein_ids)} queried."
                f"Unique taxon_ids: {len(protein_taxon_map)}"
            )

            unique_taxon_ids = list(protein_taxon_map.keys())

            if not unique_taxon_ids:
                logger.warning("No unique taxon_ids found, skipping enrichment")
                return

            # Check existing taxons
            existing_taxon_ids = await query_existing_nodes_by_ids(
                session, "Taxon", "id", unique_taxon_ids
            )
            taxon_ids_to_fetch = [tid for tid in unique_taxon_ids if tid not in existing_taxon_ids]

            logger.info(
                f"Taxon status: {len(existing_taxon_ids)} existing, {len(taxon_ids_to_fetch)} to fetch. "
                f"Total unique: {len(unique_taxon_ids)}"
            )

        logger.debug(f"Taxon IDs to fetch: {len(taxon_ids_to_fetch)}")

        # Step 2: Process existing taxons in batch
        if existing_taxon_ids:
            logger.info(f"Linking {len(existing_taxon_ids)} existing taxons to proteins")
            try:
                await self._link_existing_taxons_batch(existing_taxon_ids, protein_taxon_map)
                logger.info("Successfully linked existing taxons")
            except Exception:
                logger.exception("Failed to link existing taxons")
                raise

        # Step 3: Fetch and process new taxons in batch
        if taxon_ids_to_fetch:
            logger.info(f"Fetching {len(taxon_ids_to_fetch)} new taxons from API")
            taxonomy_responses: list[dict[str, Any]] = []
            async with httpx.AsyncClient() as client:
                async for taxonomy_response in self.taxonomy_adapter.fetch_taxa(
                    client, taxon_ids_to_fetch
                ):
                    taxonomy_responses.append(taxonomy_response)

            logger.info(f"Received {len(taxonomy_responses)} taxonomy responses from API")

            if taxonomy_responses:
                try:
                    await self._process_taxonomies_batch(taxonomy_responses, protein_taxon_map)
                    logger.info("Successfully processed new taxonomies")
                except Exception:
                    logger.exception("Failed to process new taxonomies")
                    raise

    async def _link_existing_taxons_batch(
        self, taxon_ids: list[str], protein_taxon_map: dict[str, list[str]]
    ) -> None:
        """Link existing taxons to proteins in separate transactions.

        Args:
            taxon_ids: List of existing taxon IDs to link
            protein_taxon_map: Mapping of taxon_id → [protein_ids]
        """
        if not taxon_ids:
            return

        # Collect all relationships and removals
        all_source_proteins: list[str] = []
        all_target_taxons: list[str] = []
        removal_map: dict[str, list[str]] = defaultdict(list)

        for taxon_id in taxon_ids:
            protein_ids = protein_taxon_map.get(taxon_id, [])
            if not protein_ids:
                continue
            all_source_proteins.extend(protein_ids)
            all_target_taxons.extend([taxon_id] * len(protein_ids))
            for pid in protein_ids:
                removal_map[pid].append(taxon_id)

        if not all_source_proteins:
            return

        # Transaction 1: Create relationships
        async with self._db_semaphore, self.db.async_driver.session() as session:
            await create_relationships_batch(
                session=session,
                source_label="Protein",
                source_field="id",
                source_values=all_source_proteins,
                target_label="Taxon",
                target_field="id",
                target_values=all_target_taxons,
                relationship_type="ORIGINATES_FROM",
                direction_to_source=True,
            )

        # # Transaction 2: Remove properties
        # if removal_map:
        #     async with self._db_semaphore, self.db.async_driver.session() as session:
        #         await remove_list_property_values(
        #             session=session,
        #             label="Protein",
        #             unique_field="id",
        #             unique_values=list(removal_map.keys()),
        #             list_property="taxon_ids",
        #             values_to_remove=removal_map,
        #         )

    async def _process_taxonomies_batch(
        self,
        taxonomy_responses: list[dict[str, Any]],
        protein_taxon_map: dict[str, list[str]],
    ) -> None:
        """Process multiple taxonomy responses in one transaction.

        Args:
            taxonomy_responses: List of taxonomy API responses
            protein_taxon_map: Mapping of taxon_id → [protein_ids]
        """
        if not taxonomy_responses:
            return

        # Extract data outside semaphore
        all_taxon_nodes: list[Taxon] = []
        hierarchy_infos: list[dict[str, Any]] = []
        all_source_proteins: list[str] = []
        all_target_taxons: list[str] = []
        removal_map: dict[str, list[str]] = defaultdict(list)

        for idx, taxonomy_response in enumerate(taxonomy_responses):
            try:
                # Log the response structure for debugging
                logger.debug(
                    f"Processing taxonomy response {idx + 1}/{len(taxonomy_responses)}. "
                    f"Keys: {list(taxonomy_response.keys()) if isinstance(taxonomy_response, dict) else 'not a dict'}"
                )

                hierarchy_info = self.taxonomy_adapter.extract_hierarchy_info(taxonomy_response)

                # Validate hierarchy_info structure
                if not isinstance(hierarchy_info, dict):
                    logger.error(
                        f"extract_hierarchy_info returned non-dict: {type(hierarchy_info)}. "
                        f"Response keys: {list(taxonomy_response.keys()) if isinstance(taxonomy_response, dict) else 'N/A'}"
                    )
                    continue

                # Check for required keys
                if "main_id" not in hierarchy_info:
                    logger.error(
                        f"hierarchy_info missing 'main_id' key. Keys: {list(hierarchy_info.keys())}. "
                        f"Response taxonId: {taxonomy_response.get('taxonId')}"
                    )
                    continue

                main_taxon_id = hierarchy_info["main_id"]

                if main_taxon_id is None:
                    logger.debug("Skipping taxonomy response with None main_id")
                    continue

                protein_ids = protein_taxon_map.get(str(main_taxon_id), [])
                if not protein_ids:
                    logger.debug(f"No proteins found for taxon_id {main_taxon_id}")
                    continue

                taxon_nodes = self.taxonomy_adapter.map(taxonomy_response)
                if not taxon_nodes:
                    logger.warning(f"No taxon nodes extracted for taxon_id {main_taxon_id}")
                    continue

                all_taxon_nodes.extend(taxon_nodes)
                hierarchy_infos.append(hierarchy_info)
                all_source_proteins.extend(protein_ids)
                all_target_taxons.extend([main_taxon_id] * len(protein_ids))
                for pid in protein_ids:
                    removal_map[pid].append(str(main_taxon_id))
            except KeyError:
                taxon_id = (
                    taxonomy_response.get("taxonId", "unknown")
                    if isinstance(taxonomy_response, dict)
                    else "unknown"
                )
                logger.exception(f"KeyError processing taxonomy response for taxonId {taxon_id}")
                continue
            except Exception:
                taxon_id = (
                    taxonomy_response.get("taxonId", "unknown")
                    if isinstance(taxonomy_response, dict)
                    else "unknown"
                )
                logger.exception(
                    f"Unexpected error processing taxonomy response for taxonId {taxon_id}"
                )
                continue

        if not all_taxon_nodes:
            logger.warning("No taxon nodes extracted from any taxonomy responses")
            return

        logger.info(
            f"Extracted {len(all_taxon_nodes)} taxon nodes, "
            f"{len(hierarchy_infos)} hierarchy infos, "
            f"{len(all_source_proteins)} protein relationships to create"
        )

        # Split into 4 separate transactions to reduce lock contention
        try:
            # Transaction 1: Upsert all taxon nodes
            if all_taxon_nodes:
                logger.debug(f"Upserting {len(all_taxon_nodes)} taxon nodes")
                async with self._db_semaphore, self.db.async_driver.session() as session:
                    await _upsert_nodes_with_session(session, all_taxon_nodes)

            # Transaction 2: Create all hierarchy relationships
            if hierarchy_infos:
                logger.debug(f"Creating {len(hierarchy_infos)} hierarchy chains")
                async with self._db_semaphore, self.db.async_driver.session() as session:
                    for idx, hierarchy_info in enumerate(hierarchy_infos):
                        try:
                            # Validate hierarchy_info before use
                            required_keys = ["main_id", "parent_id", "lineage_ids"]
                            missing_keys = [k for k in required_keys if k not in hierarchy_info]
                            if missing_keys:
                                logger.error(
                                    f"hierarchy_info {idx} missing keys: {missing_keys}. "
                                    f"Available keys: {list(hierarchy_info.keys())}"
                                )
                                continue

                            await create_taxonomy_hierarchy(
                                session=session,
                                main_taxon_id=hierarchy_info["main_id"],
                                parent_taxon_id=hierarchy_info["parent_id"],
                                lineage_ids=hierarchy_info["lineage_ids"],
                            )
                        except KeyError:
                            logger.exception(f"KeyError creating hierarchy {idx}")
                            raise
                        except Exception:
                            logger.exception(f"Error creating hierarchy {idx}")
                            raise

            # Transaction 3: Create all protein-taxon relationships
            if all_source_proteins:
                logger.info(f"Creating {len(all_source_proteins)} ORIGINATES_FROM relationships")
                async with self._db_semaphore, self.db.async_driver.session() as session:
                    await create_relationships_batch(
                        session=session,
                        source_label="Protein",
                        source_field="id",
                        source_values=all_source_proteins,
                        target_label="Taxon",
                        target_field="id",
                        target_values=all_target_taxons,
                        relationship_type="ORIGINATES_FROM",
                        direction_to_source=True,
                    )
                logger.info("Successfully created ORIGINATES_FROM relationships")

            # Transaction 4: Remove all taxon_ids from proteins
            if removal_map:
                logger.debug(f"Removing taxon_ids from {len(removal_map)} proteins")
                async with self._db_semaphore, self.db.async_driver.session() as session:
                    await remove_list_property_values(
                        session=session,
                        label="Protein",
                        unique_field="id",
                        unique_values=list(removal_map.keys()),
                        list_property="taxon_ids",
                        values_to_remove=removal_map,
                    )
                logger.info("Successfully removed taxon_ids from proteins")
        except KeyError:
            logger.exception(
                "KeyError in database transaction. "
                "This suggests a missing key in hierarchy_info or response data."
            )
            raise
        except Exception:
            logger.exception("Failed to persist taxonomies batch")
            raise
