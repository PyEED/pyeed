# from __future__ import annotations

# import asyncio
# from collections import defaultdict
# from typing import Any

# import httpx
# from loguru import logger
# from rich.progress import Progress, TaskID

# from pyeed.db.neo4j import GraphDB
# from pyeed.db.queries import (
#     _upsert_nodes_with_session,
#     create_reaction_molecule_relationships,
#     create_relationships_batch,
#     create_taxonomy_hierarchy,
#     query_existing_nodes_by_ids,
#     remove_list_property_values,
# )
# from pyeed.ingest.core.protocol import SENTINEL, PipelineContext
# from pyeed.ingest.model import Molecule, Reaction, Taxon
# from pyeed.ingest.sources.chebi import ChebiClient
# from pyeed.ingest.sources.rhea import RheaClient
# from pyeed.ingest.sources.taxonomy import UniProtTaxonomyAdapter


# class BatchEnrichmentStage:
#     """Batches PipelineRecords and enriches proteins with taxonomy and reactions in parallel.

#     Accumulates proteins from the pipeline, then runs targeted enrichment for taxonomy
#     and reactions concurrently without blocking the embedding pipeline.

#     Key features:
#     - Batches N proteins before enriching (default: 1000)
#     - Runs taxonomy and reaction enrichment in parallel via asyncio.gather
#     - Only enriches specific protein IDs in the batch (no full DB scan)
#     """

#     def __init__(
#         self,
#         db: GraphDB,
#         batch_size: int = 1000,
#         max_concurrent: int = 10,
#     ):
#         """Initialize batch enrichment stage.

#         Args:
#             db: GraphDB instance
#             batch_size: Number of records to accumulate before enriching
#             max_concurrent: Maximum concurrent API requests
#         """
#         self.db = db
#         self.batch_size = batch_size
#         self.max_concurrent = max_concurrent
#         self.taxonomy_adapter = UniProtTaxonomyAdapter()
#         self.rhea_client = RheaClient()
#         self.chebi_client = ChebiClient()
#         # Limit total concurrent DB operations (prevents connection pool exhaustion)
#         self._db_semaphore = asyncio.Semaphore(max_concurrent * 2)
#         self._progress_lock = asyncio.Lock()

#     async def run(
#         self,
#         input_queues: dict[str, asyncio.Queue[Any]],
#         output_queues: dict[str, asyncio.Queue[Any]],
#         context: PipelineContext,
#         progress: Progress | None,
#         task_id: TaskID | None,
#     ) -> None:
#         """Run batch enrichment: accumulate → spawn background enrichment tasks.

#         Args:
#             input_queues: Input queues with PipelineRecords
#             output_queues: Output queues (not used - terminal stage)
#             context: Pipeline context
#             progress: Progress instance
#             task_id: Task ID
#         """
#         logger.info("Starting batch enrichment stage")
#         if progress is not None and task_id is not None:
#             total = context.stats.get("total")
#             if total is not None:
#                 progress.update(task_id, total=total)

#         input_queue = next(iter(input_queues.values()))
#         batch: list[Any] = []
#         enrichment_tasks: list[asyncio.Task[Any]] = []

#         while True:
#             item = await input_queue.get()

#             # Check for total update on first item (pipeline is active by then)
#             if progress is not None and task_id is not None:
#                 total = context.stats.get("total")
#                 if total is not None:
#                     progress.update(task_id, total=total)

#             if item is SENTINEL:
#                 # Process remaining batch
#                 if batch:
#                     logger.info(f"Processing final batch of {len(batch)} records")
#                     task = asyncio.create_task(self._enrich_batch(batch.copy(), progress, task_id))
#                     enrichment_tasks.append(task)
#                     batch.clear()

#                 # Wait for all background enrichment tasks to complete
#                 if enrichment_tasks:
#                     logger.info(
#                         f"Waiting for {len(enrichment_tasks)} background enrichment tasks..."
#                     )
#                     await asyncio.gather(*enrichment_tasks, return_exceptions=True)
#                     logger.info("All enrichment tasks completed")

#                 break

#             # Accumulate records
#             batch.append(item)

#             # When batch is full, spawn background enrichment task
#             if len(batch) >= self.batch_size:
#                 logger.info(f"Batch full ({len(batch)} records), spawning enrichment task")
#                 task = asyncio.create_task(self._enrich_batch(batch.copy(), progress, task_id))
#                 enrichment_tasks.append(task)
#                 batch.clear()

#         logger.info("Batch enrichment stage complete")

#     async def _enrich_batch(
#         self,
#         batch: list[Any],
#         progress: Progress | None = None,
#         task_id: TaskID | None = None,
#     ) -> None:
#         """Enrich a batch of proteins with taxonomy and reactions in parallel.

#         Args:
#             batch: List of PipelineRecord objects
#             progress: Optional Progress instance for updating progress
#             task_id: Optional TaskID for tracking progress
#         """
#         # Extract protein IDs from batch
#         protein_ids = [record.data.id for record in batch]
#         logger.debug(f"Enriching batch of {len(protein_ids)} proteins")

#         try:
#             # Run taxonomy and reaction enrichment in parallel
#             # Each method creates its own session for concurrent execution
#             await asyncio.gather(
#                 self._enrich_taxonomy_batch(protein_ids),
#                 self._enrich_reactions_batch(protein_ids),
#                 return_exceptions=True,
#             )
#             logger.debug(f"Enrichment complete for batch of {len(protein_ids)} proteins")

#             # Advance progress now that actual work is complete
#             if progress is not None and task_id is not None:
#                 async with self._progress_lock:
#                     progress.advance(task_id, len(batch))
#         except Exception as e:
#             logger.error(f"Batch enrichment failed: {e}", exc_info=True)

#     async def _enrich_taxonomy_batch(self, protein_ids: list[str]) -> None:
#         """Enrich specific proteins with taxonomy data - streaming pattern.

#         Args:
#             protein_ids: List of protein IDs to enrich
#         """
#         # Step 1: Query proteins for taxon_ids (sequential query - single session OK)
#         async with self.db.async_driver.session() as session:
#             query = """
#             MATCH (p:Protein)
#             WHERE p.id IN $protein_ids AND p.taxon_ids IS NOT NULL AND size(p.taxon_ids) > 0
#             RETURN p.id AS id, p.taxon_ids AS taxon_ids
#             """
#             result = await session.run(query, protein_ids=protein_ids)

#             protein_taxon_map: dict[str, list[str]] = defaultdict(list)
#             async for record in result:
#                 protein_id = record["id"]
#                 taxon_ids = record["taxon_ids"]
#                 for taxon_id in taxon_ids:
#                     protein_taxon_map[taxon_id].append(protein_id)

#             unique_taxon_ids = list(protein_taxon_map.keys())

#             if not unique_taxon_ids:
#                 return

#             # Check existing taxons
#             existing_taxon_ids = await query_existing_nodes_by_ids(
#                 session, "Taxon", "id", unique_taxon_ids
#             )
#             taxon_ids_to_fetch = [tid for tid in unique_taxon_ids if tid not in existing_taxon_ids]

#         logger.debug(f"taxon ids to fetch in single batch: {taxon_ids_to_fetch}")

#         # Step 2: Process existing taxons in one batch transaction
#         if existing_taxon_ids:
#             await self._link_existing_taxons_batch(existing_taxon_ids, protein_taxon_map)

#         # Step 3: Collect new taxon API responses, then process in one batch transaction
#         if taxon_ids_to_fetch:
#             taxonomy_responses: list[dict[str, Any]] = []
#             async with httpx.AsyncClient() as client:
#                 # fetch_taxa uses batch search endpoint (no max_concurrent needed)
#                 async for taxonomy_response in self.taxonomy_adapter.fetch_taxa(
#                     client, taxon_ids_to_fetch
#                 ):
#                     taxonomy_responses.append(taxonomy_response)

#             # Process all collected responses in one transaction
#             if taxonomy_responses:
#                 await self._process_taxonomies_batch(taxonomy_responses, protein_taxon_map)

#     async def _enrich_reactions_batch(self, protein_ids: list[str]) -> None:
#         """Enrich specific proteins with reaction and molecule data - streaming pattern.

#         Args:
#             protein_ids: List of protein IDs to enrich
#         """
#         # Step 1: Query proteins for reaction_ids (sequential query - single session OK)
#         async with self.db.async_driver.session() as session:
#             query = """
#             MATCH (p:Protein)
#             WHERE p.id IN $protein_ids AND p.reaction_ids IS NOT NULL AND size(p.reaction_ids) > 0
#             RETURN p.id AS id, p.reaction_ids AS reaction_ids
#             """
#             result = await session.run(query, protein_ids=protein_ids)

#             protein_reaction_map: dict[str, list[str]] = {}
#             unique_reaction_ids: list[str] = []

#             async for record in result:
#                 protein_id = record["id"]
#                 reaction_ids = record["reaction_ids"]
#                 for reaction_id in reaction_ids:
#                     if reaction_id not in protein_reaction_map:
#                         protein_reaction_map[reaction_id] = []
#                         unique_reaction_ids.append(reaction_id)
#                     protein_reaction_map[reaction_id].append(protein_id)

#             if not unique_reaction_ids:
#                 return

#             # Check existing reactions
#             existing_reaction_ids = await query_existing_nodes_by_ids(
#                 session, "Reaction", "id", unique_reaction_ids
#             )
#             reaction_ids_to_fetch = [
#                 rid for rid in unique_reaction_ids if rid not in existing_reaction_ids
#             ]

#         # Step 2: Process existing reactions in one batch transaction
#         if existing_reaction_ids:
#             await self._link_existing_reactions_batch(existing_reaction_ids, protein_reaction_map)

#         # Step 3: Collect new reaction API responses, then process in one batch transaction
#         if reaction_ids_to_fetch:
#             reaction_data: list[tuple[dict[str, Any], dict[str, Any]]] = []
#             # fetch_reactions limits API concurrency with max_concurrent
#             async for table_row, meta_json in self.rhea_client.fetch_reactions(
#                 reaction_ids_to_fetch, max_concurrent=self.max_concurrent
#             ):
#                 reaction_data.append((table_row, meta_json))

#             # Process all collected responses in one transaction
#             if reaction_data:
#                 await self._process_reactions_batch(reaction_data, protein_reaction_map)

#     async def _process_taxonomy(
#         self,
#         taxonomy_response: dict[str, Any],
#         protein_taxon_map: dict[str, list[str]],
#     ) -> None:
#         """Process a single taxonomy response - semaphore guarded.

#         Args:
#             taxonomy_response: Taxonomy API response
#             protein_taxon_map: Mapping of taxon_id → [protein_ids]
#         """
#         # Extract data outside semaphore/lock to keep critical section small
#         try:
#             hierarchy_info = self.taxonomy_adapter.extract_hierarchy_info(taxonomy_response)
#             main_taxon_id = hierarchy_info["main_id"]

#             if main_taxon_id is None:
#                 return

#             protein_ids = protein_taxon_map.get(str(main_taxon_id), [])
#             if not protein_ids:
#                 return

#             taxon_nodes = self.taxonomy_adapter.map(taxonomy_response)
#             if not taxon_nodes:
#                 return
#         except Exception as e:
#             logger.warning(f"Failed to extract taxonomy data: {e}", exc_info=True)
#             return

#         # DB operations guarded by semaphore
#         try:
#             async with self._db_semaphore, self.db.async_driver.session() as session:
#                 await _upsert_nodes_with_session(session, taxon_nodes)

#                 await create_taxonomy_hierarchy(
#                     session=session,
#                     main_taxon_id=hierarchy_info["main_id"],
#                     parent_taxon_id=hierarchy_info["parent_id"],
#                     lineage_ids=hierarchy_info["lineage_ids"],
#                 )

#                 await create_relationships_batch(
#                     session=session,
#                     source_label="Protein",
#                     source_field="id",
#                     source_values=protein_ids,
#                     target_label="Taxon",
#                     target_field="id",
#                     target_values=[main_taxon_id] * len(protein_ids),
#                     relationship_type="ORIGINATES_FROM",
#                     direction_to_source=True,
#                 )

#                 values_to_remove = {pid: [str(main_taxon_id)] for pid in protein_ids}
#                 await remove_list_property_values(
#                     session=session,
#                     label="Protein",
#                     unique_field="id",
#                     unique_values=protein_ids,
#                     list_property="taxon_ids",
#                     values_to_remove=values_to_remove,
#                 )
#         except Exception as e:
#             taxon_id = taxonomy_response.get("taxonId")
#             logger.warning(f"Failed to persist taxonomy {taxon_id}: {e}", exc_info=True)

#     async def _process_reaction(
#         self,
#         table_row: dict[str, Any],
#         meta_json: dict[str, Any],
#         protein_reaction_map: dict[str, list[str]],
#     ) -> None:
#         """Process a single reaction response - semaphore guarded.

#         Args:
#             table_row: Rhea TSV row
#             meta_json: Rhea metadata JSON
#             protein_reaction_map: Mapping of reaction_id → [protein_ids]
#         """
#         # Extract data outside semaphore
#         try:
#             reactions = self.rhea_client._extract_reaction(table_row, meta_json, rhea_id=None)
#             if not reactions:
#                 return

#             reaction = reactions[0]  # Main reaction
#             reaction_id = reaction.id

#             protein_ids = protein_reaction_map.get(reaction_id, [])
#             # Also check for alternate IDs if primary not found directly
#             if not protein_ids:
#                 # Sometimes API returns normalized ID different from requested?
#                 # For now assume reaction.id matches one of the requested keys or we skip
#                 # Or we could try to match any ID in reactions list
#                 pass

#             # If we still don't have proteins for this reaction (maybe it was fetched
#             # as part of a batch but this specific ID isn't needed?), return
#             if not protein_ids:
#                 # It's possible the fetched reaction ID corresponds to one
#                 # requested in protein_reaction_map
#                 # But rhea_client might return multiple reactions.
#                 # We should check if ANY of the extracted reactions are in our map.
#                 pass

#             substrate_map = {reaction.id: reaction.substrate_ids}
#             product_map = {reaction.id: reaction.product_ids}
#             chebi_ids = set(reaction.substrate_ids) | set(reaction.product_ids)

#         except Exception as e:
#             logger.warning(f"Failed to extract reaction data: {e}", exc_info=True)
#             return

#         molecule_nodes: list[Molecule] = []
#         if chebi_ids:
#             # We need to check which exist. This requires DB read.
#             # We can do a small DB read here.
#             async with self._db_semaphore, self.db.async_driver.session() as session:
#                 existing_molecule_ids = await query_existing_nodes_by_ids(
#                     session, "Molecule", "id", list(chebi_ids)
#                 )

#             molecule_ids_to_fetch = [mid for mid in chebi_ids if mid not in existing_molecule_ids]

#             if molecule_ids_to_fetch:
#                 async for chebi_entry in self.chebi_client.fetch_molecules(
#                     molecule_ids_to_fetch, batch_size=50
#                 ):
#                     try:
#                         molecules = self.chebi_client._extract_molecules(chebi_entry)
#                         molecule_nodes.extend(molecules)
#                     except Exception as e:
#                         logger.warning(f"Failed to extract molecule: {e}", exc_info=True)

#         # DB operations guarded by semaphore
#         try:
#             async with self._db_semaphore, self.db.async_driver.session() as session:
#                 # Upsert nodes (Reaction + Molecules)
#                 all_nodes = [reaction, *molecule_nodes]
#                 await _upsert_nodes_with_session(session, all_nodes)

#                 # Link Reaction to Proteins
#                 if protein_ids:
#                     await create_relationships_batch(
#                         session=session,
#                         source_label="Protein",
#                         source_field="id",
#                         source_values=protein_ids,
#                         target_label="Reaction",
#                         target_field="id",
#                         target_values=[reaction.id] * len(protein_ids),
#                         relationship_type="CATALYZES",
#                         direction_to_source=True,
#                     )

#                 # Link Reaction to Molecules
#                 await create_reaction_molecule_relationships(
#                     session=session, substrate_map=substrate_map, product_map=product_map
#                 )

#                 # Remove property values
#                 if protein_ids:
#                     removal_map = {pid: [reaction.id] for pid in protein_ids}
#                     await remove_list_property_values(
#                         session=session,
#                         label="Protein",
#                         unique_field="id",
#                         unique_values=protein_ids,
#                         list_property="reaction_ids",
#                         values_to_remove=removal_map,
#                     )
#         except Exception as e:
#             logger.warning(f"Failed to persist reaction {reaction.id}: {e}", exc_info=True)

#     async def _link_existing_taxons_batch(
#         self, taxon_ids: list[str], protein_taxon_map: dict[str, list[str]]
#     ) -> None:
#         """Link multiple existing taxons to proteins in one transaction.

#         Args:
#             taxon_ids: List of taxon IDs to link
#             protein_taxon_map: Mapping of taxon_id → [protein_ids]
#         """
#         if not taxon_ids:
#             return

#         # Collect all relationships and removals
#         all_source_proteins: list[str] = []
#         all_target_taxons: list[str] = []
#         removal_map: dict[str, list[str]] = defaultdict(list)

#         for taxon_id in taxon_ids:
#             protein_ids = protein_taxon_map.get(taxon_id, [])
#             if not protein_ids:
#                 continue
#             all_source_proteins.extend(protein_ids)
#             all_target_taxons.extend([taxon_id] * len(protein_ids))
#             for pid in protein_ids:
#                 removal_map[pid].append(taxon_id)

#         if not all_source_proteins:
#             return

#         # Single transaction for all taxons
#         async with self._db_semaphore, self.db.async_driver.session() as session:
#             await create_relationships_batch(
#                 session=session,
#                 source_label="Protein",
#                 source_field="id",
#                 source_values=all_source_proteins,
#                 target_label="Taxon",
#                 target_field="id",
#                 target_values=all_target_taxons,
#                 relationship_type="ORIGINATES_FROM",
#                 direction_to_source=True,
#             )

#             await remove_list_property_values(
#                 session=session,
#                 label="Protein",
#                 unique_field="id",
#                 unique_values=list(removal_map.keys()),
#                 list_property="taxon_ids",
#                 values_to_remove=removal_map,
#             )

#     async def _process_taxonomies_batch(
#         self,
#         taxonomy_responses: list[dict[str, Any]],
#         protein_taxon_map: dict[str, list[str]],
#     ) -> None:
#         """Process multiple taxonomy responses in one transaction.

#         Args:
#             taxonomy_responses: List of taxonomy API responses
#             protein_taxon_map: Mapping of taxon_id → [protein_ids]
#         """
#         if not taxonomy_responses:
#             return

#         # Extract data outside semaphore
#         all_taxon_nodes: list[Taxon] = []
#         hierarchy_infos: list[dict[str, Any]] = []
#         all_source_proteins: list[str] = []
#         all_target_taxons: list[str] = []
#         removal_map: dict[str, list[str]] = defaultdict(list)

#         for taxonomy_response in taxonomy_responses:
#             try:
#                 hierarchy_info = self.taxonomy_adapter.extract_hierarchy_info(taxonomy_response)
#                 main_taxon_id = hierarchy_info["main_id"]

#                 if main_taxon_id is None:
#                     continue

#                 protein_ids = protein_taxon_map.get(str(main_taxon_id), [])
#                 if not protein_ids:
#                     continue

#                 taxon_nodes = self.taxonomy_adapter.map(taxonomy_response)
#                 if not taxon_nodes:
#                     continue

#                 all_taxon_nodes.extend(taxon_nodes)
#                 hierarchy_infos.append(hierarchy_info)
#                 all_source_proteins.extend(protein_ids)
#                 all_target_taxons.extend([main_taxon_id] * len(protein_ids))
#                 for pid in protein_ids:
#                     removal_map[pid].append(str(main_taxon_id))
#             except Exception as e:
#                 taxon_id = taxonomy_response.get("taxonId")
#                 logger.warning(
#                     f"Failed to extract taxonomy data for {taxon_id}: {e}", exc_info=True
#                 )

#         if not all_taxon_nodes:
#             return

#         # Single transaction for all taxonomies
#         try:
#             async with self._db_semaphore, self.db.async_driver.session() as session:
#                 # Upsert all taxon nodes
#                 await _upsert_nodes_with_session(session, all_taxon_nodes)

#                 # Create all hierarchy relationships
#                 for hierarchy_info in hierarchy_infos:
#                     await create_taxonomy_hierarchy(
#                         session=session,
#                         main_taxon_id=hierarchy_info["main_id"],
#                         parent_taxon_id=hierarchy_info["parent_id"],
#                         lineage_ids=hierarchy_info["lineage_ids"],
#                     )

#                 # Create all protein-taxon relationships
#                 if all_source_proteins:
#                     await create_relationships_batch(
#                         session=session,
#                         source_label="Protein",
#                         source_field="id",
#                         source_values=all_source_proteins,
#                         target_label="Taxon",
#                         target_field="id",
#                         target_values=all_target_taxons,
#                         relationship_type="ORIGINATES_FROM",
#                         direction_to_source=True,
#                     )

#                 # Remove all taxon_ids from proteins
#                 if removal_map:
#                     await remove_list_property_values(
#                         session=session,
#                         label="Protein",
#                         unique_field="id",
#                         unique_values=list(removal_map.keys()),
#                         list_property="taxon_ids",
#                         values_to_remove=removal_map,
#                     )
#         except Exception as e:
#             logger.warning(f"Failed to persist taxonomies batch: {e}", exc_info=True)

#     async def _link_existing_reactions_batch(
#         self, reaction_ids: list[str], protein_reaction_map: dict[str, list[str]]
#     ) -> None:
#         """Link multiple existing reactions to proteins and enrich with molecules.

#         Args:
#             reaction_ids: List of reaction IDs to link
#             protein_reaction_map: Mapping of reaction_id → [protein_ids]
#         """
#         if not reaction_ids:
#             return

#         # Collect all relationships and removals
#         all_source_proteins: list[str] = []
#         all_target_reactions: list[str] = []
#         removal_map: dict[str, list[str]] = defaultdict(list)

#         for reaction_id in reaction_ids:
#             protein_ids = protein_reaction_map.get(reaction_id, [])
#             if not protein_ids:
#                 continue
#             all_source_proteins.extend(protein_ids)
#             all_target_reactions.extend([reaction_id] * len(protein_ids))
#             for pid in protein_ids:
#                 removal_map[pid].append(reaction_id)

#         if not all_source_proteins:
#             return

#         # Query existing Reaction nodes for their substrate_ids and product_ids
#         async with self.db.async_driver.session() as session:
#             query = """
#             MATCH (r:Reaction)
#             WHERE r.id IN $reaction_ids
#             RETURN r.id AS id, r.substrate_ids AS substrate_ids, r.product_ids AS product_ids
#             """
#             result = await session.run(query, reaction_ids=reaction_ids)

#             substrate_map: dict[str, list[str]] = {}
#             product_map: dict[str, list[str]] = {}
#             all_chebi_ids: set[str] = set()

#             async for record in result:
#                 reaction_id = record["id"]
#                 substrate_ids = record.get("substrate_ids") or []
#                 product_ids = record.get("product_ids") or []

#                 if substrate_ids:
#                     substrate_map[reaction_id] = substrate_ids
#                     all_chebi_ids.update(substrate_ids)
#                 if product_ids:
#                     product_map[reaction_id] = product_ids
#                     all_chebi_ids.update(product_ids)

#         # Fetch molecules for all reactions
#         all_molecule_nodes: list[Molecule] = []
#         if all_chebi_ids:
#             async with self._db_semaphore, self.db.async_driver.session() as session:
#                 existing_molecule_ids = await query_existing_nodes_by_ids(
#                     session, "Molecule", "id", list(all_chebi_ids)
#                 )

#             molecule_ids_to_fetch = [
#                 mid for mid in all_chebi_ids if mid not in existing_molecule_ids
#             ]

#             if molecule_ids_to_fetch:
#                 async for chebi_entry in self.chebi_client.fetch_molecules(
#                     molecule_ids_to_fetch, batch_size=50
#                 ):
#                     try:
#                         molecules = self.chebi_client._extract_molecules(chebi_entry)
#                         all_molecule_nodes.extend(molecules)
#                     except Exception as e:
#                         logger.warning(f"Failed to extract molecule: {e}", exc_info=True)

#         # Single transaction for all operations
#         try:
#             async with self._db_semaphore, self.db.async_driver.session() as session:
#                 # Upsert molecules if any were fetched
#                 if all_molecule_nodes:
#                     await _upsert_nodes_with_session(session, all_molecule_nodes)

#                 # Link Reactions to Proteins
#                 await create_relationships_batch(
#                     session=session,
#                     source_label="Protein",
#                     source_field="id",
#                     source_values=all_source_proteins,
#                     target_label="Reaction",
#                     target_field="id",
#                     target_values=all_target_reactions,
#                     relationship_type="CATALYZES",
#                     direction_to_source=True,
#                     tx_size=5000,
#                 )

#                 # Link Reactions to Molecules
#                 if substrate_map or product_map:
#                     await create_reaction_molecule_relationships(
#                         session=session, substrate_map=substrate_map, product_map=product_map
#                     )

#                 # Remove reaction_ids from proteins
#                 await remove_list_property_values(
#                     session=session,
#                     label="Protein",
#                     unique_field="id",
#                     unique_values=list(removal_map.keys()),
#                     list_property="reaction_ids",
#                     values_to_remove=removal_map,
#                 )
#         except Exception as e:
#             logger.warning(f"Failed to persist existing reactions batch: {e}", exc_info=True)

#     async def _process_reactions_batch(
#         self,
#         reaction_data: list[tuple[dict[str, Any], dict[str, Any]]],
#         protein_reaction_map: dict[str, list[str]],
#     ) -> None:
#         """Process multiple reaction responses in one transaction.

#         Args:
#             reaction_data: List of (table_row, meta_json) tuples from Rhea API
#             protein_reaction_map: Mapping of reaction_id → [protein_ids]
#         """
#         if not reaction_data:
#             return

#         # Extract data outside semaphore
#         all_reaction_nodes: list[Reaction] = []
#         all_molecule_nodes: list[Molecule] = []
#         substrate_map: dict[str, list[str]] = {}
#         product_map: dict[str, list[str]] = {}
#         all_source_proteins: list[str] = []
#         all_target_reactions: list[str] = []
#         removal_map: dict[str, list[str]] = defaultdict(list)
#         all_chebi_ids: set[str] = set()

#         for table_row, meta_json in reaction_data:
#             try:
#                 reactions = self.rhea_client._extract_reaction(table_row, meta_json, rhea_id=None)
#                 if not reactions:
#                     continue

#                 reaction = reactions[0]  # Main reaction
#                 reaction_id = reaction.id

#                 protein_ids = protein_reaction_map.get(reaction_id, [])
#                 if not protein_ids:
#                     continue

#                 all_reaction_nodes.append(reaction)
#                 substrate_map[reaction.id] = reaction.substrate_ids
#                 product_map[reaction.id] = reaction.product_ids
#                 all_chebi_ids.update(reaction.substrate_ids)
#                 all_chebi_ids.update(reaction.product_ids)
#                 all_source_proteins.extend(protein_ids)
#                 all_target_reactions.extend([reaction.id] * len(protein_ids))
#                 for pid in protein_ids:
#                     removal_map[pid].append(reaction.id)
#             except Exception as e:
#                 logger.warning(f"Failed to extract reaction data: {e}", exc_info=True)

#         if not all_reaction_nodes:
#             return

#         # Fetch molecules for all reactions
#         if all_chebi_ids:
#             async with self._db_semaphore, self.db.async_driver.session() as session:
#                 existing_molecule_ids = await query_existing_nodes_by_ids(
#                     session, "Molecule", "id", list(all_chebi_ids)
#                 )

#             molecule_ids_to_fetch = [
#                 mid for mid in all_chebi_ids if mid not in existing_molecule_ids
#             ]

#             if molecule_ids_to_fetch:
#                 async for chebi_entry in self.chebi_client.fetch_molecules(
#                     molecule_ids_to_fetch, batch_size=50
#                 ):
#                     try:
#                         molecules = self.chebi_client._extract_molecules(chebi_entry)
#                         all_molecule_nodes.extend(molecules)
#                     except Exception as e:
#                         logger.warning(f"Failed to extract molecule: {e}", exc_info=True)

#         # Single transaction for all reactions and molecules
#         try:
#             async with self._db_semaphore, self.db.async_driver.session() as session:
#                 # Upsert all nodes (Reactions + Molecules)
#                 all_nodes = all_reaction_nodes + all_molecule_nodes
#                 if all_nodes:
#                     await _upsert_nodes_with_session(session, all_nodes)

#                 # Link all Reactions to Proteins
#                 if all_source_proteins:
#                     await create_relationships_batch(
#                         session=session,
#                         source_label="Protein",
#                         source_field="id",
#                         source_values=all_source_proteins,
#                         target_label="Reaction",
#                         target_field="id",
#                         target_values=all_target_reactions,
#                         relationship_type="CATALYZES",
#                         direction_to_source=True,
#                     )

#                 # Link all Reactions to Molecules
#                 if substrate_map or product_map:
#                     await create_reaction_molecule_relationships(
#                         session=session, substrate_map=substrate_map, product_map=product_map
#                     )

#                 # Remove all reaction_ids from proteins
#                 if removal_map:
#                     await remove_list_property_values(
#                         session=session,
#                         label="Protein",
#                         unique_field="id",
#                         unique_values=list(removal_map.keys()),
#                         list_property="reaction_ids",
#                         values_to_remove=removal_map,
#                     )
#         except Exception as e:
#             logger.warning(f"Failed to persist reactions batch: {e}", exc_info=True)
