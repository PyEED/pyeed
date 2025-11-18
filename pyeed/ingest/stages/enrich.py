"""Enrichment stages for multi-stage pipeline processing."""

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
    create_reaction_molecule_relationships,
    create_relationships_batch,
    create_taxonomy_hierarchy,
    query_existing_nodes_by_ids,
    query_nodes_by_list_property,
    remove_list_property_values,
)
from pyeed.ingest.core.protocol import PipelineContext
from pyeed.ingest.model.molecule import Molecule
from pyeed.ingest.model.reaction import Reaction
from pyeed.ingest.model.taxon import Taxon
from pyeed.ingest.sources.chebi import ChebiClient
from pyeed.ingest.sources.rhea import RheaClient
from pyeed.ingest.sources.taxonomy import UniProtTaxonomyAdapter


class TaxonomyEnrichmentStage:
    """Enriches Protein nodes with taxonomy lineage.

    Queries DB for Proteins with taxon_ids, fetches taxonomy data concurrently,
    creates Taxon nodes with IS_A hierarchy, links via ORIGINATES_FROM.

    Example:
        For a Protein with taxon_ids=["9606"]:
        1. Fetches taxonomy for 9606 (Homo sapiens)
        2. Upserts Taxon nodes (main + parent + lineage)
        3. Creates IS_A hierarchy: Homo sapiens IS_A Homo IS_A ... IS_A cellular organisms
        4. Creates ORIGINATES_FROM: Protein → Homo sapiens
        5. Removes "9606" from protein.taxon_ids
    """

    def __init__(
        self,
        db: GraphDB,
        max_concurrent: int = 10,
        batch_size: int = 50,
        ensure_schema: bool = True,
    ):
        """Initialize taxonomy enrichment stage.

        Args:
            db: GraphDB instance
            max_concurrent: Maximum concurrent API requests
            batch_size: Batch size for processing (currently unused)
            ensure_schema: If True, sync schema constraints before enrichment
        """
        self.db = db
        self.max_concurrent = max_concurrent
        self.batch_size = batch_size
        self.ensure_schema = ensure_schema
        self.adapter = UniProtTaxonomyAdapter()

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Run enrichment: query → fetch → upsert → link → cleanup.

        This stage queries the database directly and updates it, so it doesn't
        consume from input queues or produce to output queues.

        Args:
            input_queues: Named input queues (empty for this stage)
            output_queues: Named output queues (empty for this stage)
            context: Shared pipeline context
            progress: Progress object for reporting
            task_id: TaskID for progress tracking
        """
        logger.info("Starting taxonomy enrichment stage")

        # Ensure schema constraints exist for Taxon nodes
        if self.ensure_schema:
            await self.db.sync_schema([Taxon])
            logger.debug("Schema constraints synced for Taxon")

        async with self.db.async_driver.session() as session:
            # Step 1: Query Proteins with taxon_ids
            logger.info("Querying Proteins with taxon_ids")
            protein_taxon_map: dict[str, list[str]] = defaultdict(list)

            async for record in query_nodes_by_list_property(
                session=session,
                label="Protein",
                list_property="taxon_ids",
                unique_field="id",
            ):
                protein_id = record["id"]
                taxon_ids = record["taxon_ids"]

                # Build mapping: taxon_id → [protein_ids]
                for taxon_id in taxon_ids:
                    protein_taxon_map[taxon_id].append(protein_id)

            unique_taxon_ids = list(protein_taxon_map.keys())
            logger.info(
                f"Found {len(protein_taxon_map)} unique taxon_ids "
                f"across proteins needing enrichment"
            )

            if not unique_taxon_ids:
                logger.info("No proteins need taxonomy enrichment")
                return

            # Step 2: Check which taxons already exist in database
            existing_taxon_ids = await self._get_existing_taxon_ids(session, unique_taxon_ids)
            taxon_ids_to_fetch = [tid for tid in unique_taxon_ids if tid not in existing_taxon_ids]

            logger.info(
                f"Found {len(existing_taxon_ids)} existing taxa, "
                f"fetching {len(taxon_ids_to_fetch)} new taxa"
            )

            # Step 3: Link proteins to existing taxons (no API fetch needed)
            for taxon_id in existing_taxon_ids:
                await self._link_existing_taxon(session, taxon_id, protein_taxon_map)
                if progress and task_id:
                    progress.update(task_id, advance=1)

            # Step 4: Fetch and process new taxons from API
            enriched_count = len(existing_taxon_ids)
            failed_count = 0

            if taxon_ids_to_fetch:
                async with httpx.AsyncClient() as client:
                    async for taxonomy_response in self.adapter.fetch_taxa(
                        client, taxon_ids_to_fetch, max_concurrent=self.max_concurrent
                    ):
                        try:
                            await self._process_taxonomy(
                                session=session,
                                taxonomy_response=taxonomy_response,
                                protein_taxon_map=protein_taxon_map,
                            )
                            enriched_count += 1

                            if progress and task_id:
                                progress.update(task_id, advance=1)

                        except Exception as e:
                            taxon_id = taxonomy_response.get("taxonId")
                            logger.warning(
                                f"Failed to process taxonomy {taxon_id}: {e}",
                                exc_info=True,
                            )
                            failed_count += 1

        logger.info(
            f"Taxonomy enrichment complete: {enriched_count} enriched, {failed_count} failed"
        )


class ReactionEnrichmentStage:
    """Enriches Protein nodes with reaction and molecule information.

    Queries DB for Proteins with reaction_ids, fetches reaction data from Rhea,
    fetches molecule data from ChEBI, creates Reaction and Molecule nodes,
    links via CATALYZES, HAS_SUBSTRATE, HAS_PRODUCT relationships.

    Example:
        For a Protein with reaction_ids=["RHEA:13065"]:
        1. Fetches reaction from Rhea API (13065)
        2. Extracts substrate/product ChEBI IDs from reaction
        3. Fetches molecule data from ChEBI API
        4. Upserts Reaction and Molecule nodes
        5. Creates CATALYZES: Protein → Reaction
        6. Creates HAS_SUBSTRATE: Reaction → Molecule (substrates)
        7. Creates HAS_PRODUCT: Reaction → Molecule (products)
        8. Removes "RHEA:13065" from protein.reaction_ids
    """

    def __init__(
        self,
        db: GraphDB,
        max_concurrent: int = 10,
        batch_size: int = 50,
        ensure_schema: bool = True,
    ):
        """Initialize reaction enrichment stage.

        Args:
            db: Neo4j database connection
            max_concurrent: Maximum concurrent API requests
            batch_size: Batch size for processing
            ensure_schema: Whether to sync schema at startup
        """
        self.db = db
        self.max_concurrent = max_concurrent
        self.batch_size = batch_size
        self.ensure_schema = ensure_schema
        self.rhea_client = RheaClient()
        self.chebi_client = ChebiClient()

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Run reaction enrichment stage.

        Steps:
        1. Query proteins with reaction_ids from database
        2. Check which reactions already exist
        3. Link proteins to existing reactions
        4. Fetch new reactions from Rhea API
        5. Extract ChEBI IDs from reactions
        6. Check which molecules already exist
        7. Fetch new molecules from ChEBI API
        8. Upsert Reaction and Molecule nodes
        9. Create relationships (CATALYZES, HAS_SUBSTRATE, HAS_PRODUCT)
        10. Remove processed reaction_ids from proteins
        """
        async with self.db.async_driver.session() as session:
            # Ensure schema constraints exist for Reaction and Molecule nodes
            if self.ensure_schema:
                await self.db.sync_schema([Reaction, Molecule])
                logger.debug("Schema constraints synced for Reaction and Molecule")

            # Step 1: Query Proteins with reaction_ids
            logger.info("Querying proteins with reaction_ids...")
            protein_reaction_map: dict[str, list[str]] = {}
            unique_reaction_ids: list[str] = []

            async for record in query_nodes_by_list_property(
                session=session,
                label="Protein",
                list_property="reaction_ids",
                unique_field="id",
            ):
                protein_id = record["id"]
                reaction_ids = record["reaction_ids"]

                for reaction_id in reaction_ids:
                    if reaction_id not in protein_reaction_map:
                        protein_reaction_map[reaction_id] = []
                        unique_reaction_ids.append(reaction_id)
                    protein_reaction_map[reaction_id].append(protein_id)

            if not unique_reaction_ids:
                logger.info("No proteins with reaction_ids found")
                return

            logger.info(
                f"Found {len(unique_reaction_ids)} unique reaction IDs across "
                f"{sum(len(v) for v in protein_reaction_map.values())} protein-reaction links"
            )

            # Update progress total
            if progress and task_id:
                progress.update(task_id, total=len(unique_reaction_ids))

            # Step 2: Check existing Reaction nodes
            existing_reaction_ids = await query_existing_nodes_by_ids(
                session, "Reaction", "id", unique_reaction_ids
            )
            reaction_ids_to_fetch = [
                rid for rid in unique_reaction_ids if rid not in existing_reaction_ids
            ]

            logger.info(
                f"{len(existing_reaction_ids)} reactions already exist, "
                f"{len(reaction_ids_to_fetch)} need to be fetched"
            )

            # Step 3: Link proteins to existing reactions
            for reaction_id in existing_reaction_ids:
                await self._link_existing_reaction(session, reaction_id, protein_reaction_map)
                if progress and task_id:
                    progress.update(task_id, advance=1)

            # Step 4: Fetch new reactions from Rhea API
            if not reaction_ids_to_fetch:
                logger.info("All reactions already exist, skipping API fetch")
                return

            logger.info(f"Fetching {len(reaction_ids_to_fetch)} reactions from Rhea API...")

            # Collect reaction data and extract ChEBI IDs
            reaction_nodes: list[Reaction] = []
            substrate_map: dict[str, list[str]] = {}
            product_map: dict[str, list[str]] = {}
            all_chebi_ids: set[str] = set()
            enriched_count = 0
            failed_count = 0

            async with httpx.AsyncClient() as http_client:
                async for table_row, meta_json in self.rhea_client.fetch_reactions(
                    reaction_ids_to_fetch, max_concurrent=self.max_concurrent
                ):
                    try:
                        # Extract reaction using existing map method
                        reactions = self.rhea_client._extract_reaction(
                            table_row, meta_json, rhea_id=None
                        )

                        if not reactions:
                            logger.debug(f"No reaction data extracted from {table_row}")
                            failed_count += 1
                            continue

                        reaction = reactions[0]  # Should only be one
                        reaction_nodes.append(reaction)

                        # Store substrate/product maps
                        substrate_map[reaction.id] = reaction.substrate_ids
                        product_map[reaction.id] = reaction.product_ids

                        # Collect all ChEBI IDs
                        all_chebi_ids.update(reaction.substrate_ids)
                        all_chebi_ids.update(reaction.product_ids)

                        enriched_count += 1

                    except Exception as e:
                        logger.warning(f"Failed to process reaction: {e}", exc_info=True)
                        failed_count += 1

            logger.info(
                f"Extracted {len(reaction_nodes)} reactions with "
                f"{len(all_chebi_ids)} unique molecule IDs"
            )

            # Step 5-6: Check existing Molecule nodes and fetch new ones
            if all_chebi_ids:
                chebi_id_list = list(all_chebi_ids)
                existing_molecule_ids = await query_existing_nodes_by_ids(
                    session, "Molecule", "id", chebi_id_list
                )
                molecule_ids_to_fetch = [
                    mid for mid in chebi_id_list if mid not in existing_molecule_ids
                ]

                logger.info(
                    f"{len(existing_molecule_ids)} molecules already exist, "
                    f"{len(molecule_ids_to_fetch)} need to be fetched"
                )

                # Step 7: Fetch new molecules from ChEBI API
                molecule_nodes: list[Molecule] = []
                if molecule_ids_to_fetch:
                    logger.info(
                        f"Fetching {len(molecule_ids_to_fetch)} molecules from ChEBI API..."
                    )

                    async for chebi_entry in self.chebi_client.fetch_molecules(
                        molecule_ids_to_fetch, batch_size=self.batch_size
                    ):
                        try:
                            molecules = self.chebi_client._extract_molecules(chebi_entry)
                            molecule_nodes.extend(molecules)
                        except Exception as e:
                            logger.warning(f"Failed to extract molecule: {e}", exc_info=True)

                    logger.info(f"Fetched {len(molecule_nodes)} molecules from ChEBI API")

            else:
                molecule_nodes = []

            # Step 8: Upsert Reaction and Molecule nodes
            all_nodes = reaction_nodes + molecule_nodes
            if all_nodes:
                logger.info(
                    f"Upserting {len(reaction_nodes)} reactions and {len(molecule_nodes)} molecules..."
                )
                await _upsert_nodes_with_session(session, all_nodes)

            # Step 9: Create relationships
            # Protein -CATALYZES-> Reaction
            for reaction_id, protein_ids in protein_reaction_map.items():
                if (
                    reaction_id in [r.id for r in reaction_nodes]
                    or reaction_id in existing_reaction_ids
                ):
                    await create_relationships_batch(
                        session=session,
                        source_label="Protein",
                        source_field="id",
                        source_values=protein_ids,
                        target_label="Reaction",
                        target_field="id",
                        target_values=[reaction_id] * len(protein_ids),
                        relationship_type="CATALYZES",
                        direction_to_source=True,
                        tx_size=5000,
                    )

            # Reaction -HAS_SUBSTRATE/PRODUCT-> Molecule
            await create_reaction_molecule_relationships(
                session=session,
                substrate_map=substrate_map,
                product_map=product_map,
            )

            # Step 10: Remove processed reaction_ids from proteins
            removal_map: dict[str, list[str]] = {}
            for reaction_id, protein_ids in protein_reaction_map.items():
                for protein_id in protein_ids:
                    if protein_id not in removal_map:
                        removal_map[protein_id] = []
                    removal_map[protein_id].append(reaction_id)

            await remove_list_property_values(
                session=session,
                label="Protein",
                unique_field="id",
                unique_values=list(removal_map.keys()),
                list_property="reaction_ids",
                values_to_remove=removal_map,
            )

            if progress and task_id:
                progress.update(task_id, advance=len(reaction_ids_to_fetch))

        logger.info(
            f"Reaction enrichment complete: {enriched_count} enriched, {failed_count} failed"
        )

    async def _link_existing_reaction(
        self, session: Any, reaction_id: str, protein_reaction_map: dict[str, list[str]]
    ) -> None:
        """Link proteins to an existing reaction and remove from protein.reaction_ids.

        Args:
            session: Neo4j session
            reaction_id: Reaction ID that already exists
            protein_reaction_map: Mapping of reaction_id -> [protein_ids]
        """
        protein_ids = protein_reaction_map.get(reaction_id, [])
        if not protein_ids:
            return

        # Create Protein -CATALYZES-> Reaction relationships
        await create_relationships_batch(
            session=session,
            source_label="Protein",
            source_field="id",
            source_values=protein_ids,
            target_label="Reaction",
            target_field="id",
            target_values=[reaction_id] * len(protein_ids),
            relationship_type="CATALYZES",
            direction_to_source=True,
            tx_size=5000,
        )

        # Remove reaction_id from protein.reaction_ids
        removal_map = {protein_id: [reaction_id] for protein_id in protein_ids}
        await remove_list_property_values(
            session=session,
            label="Protein",
            unique_field="id",
            unique_values=protein_ids,
            list_property="reaction_ids",
            values_to_remove=removal_map,
        )

    async def _process_taxonomy(
        self,
        session: Any,
        taxonomy_response: dict[str, Any],
        protein_taxon_map: dict[str, list[str]],
    ) -> None:
        """Process a single taxonomy response: upsert, link, cleanup.

        Args:
            session: Neo4j async session
            taxonomy_response: Taxonomy API response
            protein_taxon_map: Mapping of taxon_id → [protein_ids]
        """
        # Extract hierarchy info
        hierarchy_info = self.adapter.extract_hierarchy_info(taxonomy_response)
        main_taxon_id = hierarchy_info["main_id"]

        if main_taxon_id is None:
            logger.warning("Taxonomy response has no main taxon ID, skipping")
            return

        # Get proteins that need this taxon
        protein_ids = protein_taxon_map.get(str(main_taxon_id), [])
        if not protein_ids:
            logger.debug(f"No proteins waiting for taxon {main_taxon_id}, skipping")
            return

        # Parse into Taxon objects
        taxon_nodes = self.adapter.map(taxonomy_response)
        if not taxon_nodes:
            logger.warning(f"Failed to parse taxonomy {main_taxon_id}, skipping")
            return

        logger.debug(
            f"Processing taxon {main_taxon_id} with {len(taxon_nodes)} "
            f"nodes for {len(protein_ids)} proteins"
        )

        # Step 1: Upsert all Taxon nodes
        await _upsert_nodes_with_session(session, taxon_nodes)

        # Step 2: Create IS_A hierarchy
        await create_taxonomy_hierarchy(
            session=session,
            main_taxon_id=hierarchy_info["main_id"],
            parent_taxon_id=hierarchy_info["parent_id"],
            lineage_ids=hierarchy_info["lineage_ids"],
        )

        # Step 3: Create ORIGINATES_FROM relationships (Protein → main taxon)
        await create_relationships_batch(
            session=session,
            source_label="Protein",
            source_field="id",
            source_values=protein_ids,
            target_label="Taxon",
            target_field="id",
            target_values=[main_taxon_id] * len(protein_ids),
            relationship_type="ORIGINATES_FROM",
            direction_to_source=True,  # Protein -> Taxon
        )

        # Step 4: Remove taxon_id from protein.taxon_ids
        values_to_remove = {pid: [str(main_taxon_id)] for pid in protein_ids}
        await remove_list_property_values(
            session=session,
            label="Protein",
            unique_field="id",
            unique_values=protein_ids,
            list_property="taxon_ids",
            values_to_remove=values_to_remove,
        )

        logger.debug(f"Enriched {len(protein_ids)} proteins with taxon {main_taxon_id}")

    async def _get_existing_taxon_ids(self, session: Any, taxon_ids: list[str]) -> set[str]:
        """Check which taxon IDs already exist in database.

        Args:
            session: Neo4j async session
            taxon_ids: List of taxon IDs (as strings) to check

        Returns:
            Set of taxon IDs (as strings) that exist in database
        """
        if not taxon_ids:
            return set()

        query = """
        MATCH (t:Taxon)
        WHERE t.id IN $taxon_ids
        RETURN t.id AS taxon_id
        """
        result = await session.run(query, taxon_ids=taxon_ids)
        existing = {str(record["taxon_id"]) async for record in result}
        return existing

    async def _link_existing_taxon(
        self, session: Any, taxon_id: str, protein_taxon_map: dict[str, list[str]]
    ) -> None:
        """Link proteins to existing taxon without fetching from API.

        Args:
            session: Neo4j async session
            taxon_id: Taxon ID (as string)
            protein_taxon_map: Mapping of taxon_id → [protein_ids]
        """
        protein_ids = protein_taxon_map.get(taxon_id, [])
        if not protein_ids:
            return

        # Create ORIGINATES_FROM relationships
        await create_relationships_batch(
            session=session,
            source_label="Protein",
            source_field="id",
            source_values=protein_ids,
            target_label="Taxon",
            target_field="id",
            target_values=[taxon_id] * len(protein_ids),
            relationship_type="ORIGINATES_FROM",
            direction_to_source=True,
        )

        # Remove taxon_id from protein.taxon_ids
        values_to_remove = {pid: [taxon_id] for pid in protein_ids}
        await remove_list_property_values(
            session=session,
            label="Protein",
            unique_field="id",
            unique_values=protein_ids,
            list_property="taxon_ids",
            values_to_remove=values_to_remove,
        )

        logger.debug(f"Linked {len(protein_ids)} proteins to existing taxon {taxon_id}")
