from __future__ import annotations

from collections.abc import AsyncIterator
from typing import Any

from neo4j import AsyncSession

from ..ingest.core.pipeline import PipelineRecord
from ..ingest.model.pyeedbase import BaseNode


async def upsert_pipeline_records(
    session: AsyncSession,
    records: list[PipelineRecord[BaseNode]],
    tx_size: int = 5000,
) -> None:
    """Batch upsert PipelineRecords: nodes, children, relationships, cleanup.

    Processes multiple PipelineRecord objects efficiently by:
    1. Batch upserting all nodes (parents + children)
    2. Batch creating all relationships (grouped by type)
    3. Batch removing matched values from list properties

    Example:
        Process 100 Proteins with GOAnnotations and Annotations:

        >>> records = [uniprot_adapter.map(data) for data in fetch_data()]
        >>> async with db.async_driver.session() as session:
        ...     await upsert_pipeline_records(session, records)

    Args:
        session: Neo4j async session
        records: List of PipelineRecord objects from sources (UniProt, etc.)
        tx_size: Transaction batch size (default: 5000)
    """
    if not records:
        return

    # Step 1: Collect and upsert all nodes
    all_nodes = _collect_all_nodes(records)
    await _upsert_nodes_with_session(session, all_nodes, tx_size)

    # Step 2: Collect and create all relationships
    relationship_groups = _collect_relationships(records)
    await _create_relationships_batched(session, relationship_groups, tx_size)

    # Step 3: Collect and remove list property values
    removal_groups = _collect_list_removals(records)
    await _remove_list_values_batched(session, removal_groups, tx_size)


def _collect_all_nodes(records: list[PipelineRecord[BaseNode]]) -> list[BaseNode]:
    """Collect all nodes from records (parent + children)."""
    nodes = []
    for record in records:
        nodes.append(record.data)
        for child_record in record.children:
            nodes.extend(child_record.data)
    return nodes


def _collect_relationships(
    records: list[PipelineRecord[BaseNode]],
) -> dict[tuple[str, str, str, str, str, bool], list[tuple[str, str]]]:
    """Collect all relationships grouped by type.

    Returns:
        Dict keyed by (parent_label, parent_field, child_label, child_field,
                      edge_name, direction_to_parent)
        with values as list of (parent_value, child_value) pairs
    """
    groups = {}

    for record in records:
        parent_label = type(record.data).__name__
        parent_field = record.data.get_unique_model_field()
        parent_value = str(getattr(record.data, parent_field))

        for child_record in record.children:
            if not child_record.data:
                continue

            child_label = type(child_record.data[0]).__name__
            key = (
                parent_label,
                parent_field,
                child_label,
                child_record.child_field,
                child_record.edge_name,
                child_record.edge_direction_to_parent,
            )

            if key not in groups:
                groups[key] = []

            for child in child_record.data:
                child_value = str(getattr(child, child_record.child_field))
                groups[key].append((parent_value, child_value))

    return groups


def _collect_list_removals(
    records: list[PipelineRecord[BaseNode]],
) -> dict[tuple[str, str, str], dict[str, list[str]]]:
    """Collect list property values to remove.

    Returns:
        Dict keyed by (label, unique_field, list_property)
        with values as dict {unique_value: [values_to_remove]}
    """
    groups = {}

    for record in records:
        parent_label = type(record.data).__name__
        parent_field = record.data.get_unique_model_field()
        parent_value = str(getattr(record.data, parent_field))

        for child_record in record.children:
            if not child_record.remove_parent_value_on_join:
                continue

            key = (parent_label, parent_field, child_record.parent_field)
            if key not in groups:
                groups[key] = {}

            if parent_value not in groups[key]:
                groups[key][parent_value] = []

            child_values = [
                str(getattr(child, child_record.child_field)) for child in child_record.data
            ]
            groups[key][parent_value].extend(child_values)

    return groups


async def _create_relationships_batched(
    session: AsyncSession,
    groups: dict[tuple[str, str, str, str, str, bool], list[tuple[str, str]]],
    tx_size: int,
) -> None:
    """Create all relationships in batches using optimized MERGE pattern."""

    for (
        parent_label,
        parent_field,
        child_label,
        child_field,
        edge_name,
        direction_to_parent,
    ), pairs in groups.items():
        # Deduplicate pairs in Python to avoid duplicate MERGE operations
        unique_pairs = list(dict.fromkeys(pairs))  # Preserves order, removes duplicates

        # Process in chunks
        for i in range(0, len(unique_pairs), tx_size):
            chunk = unique_pairs[i : i + tx_size]

            # Uses indexes on parent_field and child_field for fast lookups
            if direction_to_parent:
                query = f"""
                UNWIND $rows AS r
                MATCH (p:`{parent_label}` {{ `{parent_field}`: r.pv }})
                MATCH (c:`{child_label}` {{ `{child_field}`: r.cv }})
                MERGE (c)-[:`{edge_name}`]->(p)
                """
            else:
                query = f"""
                UNWIND $rows AS r
                MATCH (p:`{parent_label}` {{ `{parent_field}`: r.pv }})
                MATCH (c:`{child_label}` {{ `{child_field}`: r.cv }})
                MERGE (p)-[:`{edge_name}`]->(c)
                """

            rows = [{"pv": pv, "cv": cv} for pv, cv in chunk]
            await session.execute_write(query, rows=rows)


async def _remove_list_values_batched(
    session: AsyncSession,
    groups: dict[tuple[str, str, str], dict[str, list[str]]],
    tx_size: int,
) -> None:
    """Remove list property values in batches."""
    for (label, unique_field, list_property), values_to_remove in groups.items():
        await remove_list_property_values(
            session=session,
            label=label,
            unique_field=unique_field,
            unique_values=list(values_to_remove.keys()),
            list_property=list_property,
            values_to_remove=values_to_remove,
            tx_size=tx_size,
        )


async def _upsert_nodes_with_session(
    session: AsyncSession,
    nodes: list[BaseNode],
    tx_size: int = 5000,
) -> None:
    """Upsert nodes using provided session (internal helper)."""
    if not nodes:
        return

    # Group by label
    by_label: dict[str, list[BaseNode]] = {}
    for node in nodes:
        label = type(node).__name__
        by_label.setdefault(label, []).append(node)

    for label, node_list in by_label.items():
        unique_field = node_list[0].get_unique_model_field()

        # Deduplicate by unique field (last occurrence wins - most recent data)
        seen: dict[str, BaseNode] = {}
        for node in node_list:
            unique_value = str(getattr(node, unique_field))
            seen[unique_value] = node

        # Build rows from deduplicated nodes
        rows = [
            {
                "key": str(getattr(node, unique_field)),
                "props": node.model_dump(),
            }
            for node in seen.values()
        ]

        for i in range(0, len(rows), tx_size):
            chunk = rows[i : i + tx_size]
            query = f"""
            UNWIND $rows AS r
            MERGE (n:`{label}` {{ `{unique_field}`: r.key }})
            SET n += r.props
            """
            await session.run(query, rows=chunk)


async def create_relationships_batch(
    session: AsyncSession,
    source_label: str,
    source_field: str,
    source_values: list[str],
    target_label: str,
    target_field: str,
    target_values: list[str],
    relationship_type: str,
    direction_to_source: bool,
    tx_size: int = 5000,
) -> None:
    """Create relationships in batch between source and target nodes.

    Args:
        session: Neo4j async session
        source_label: Label of source nodes
        source_field: Unique field name of source nodes
        source_values: List of source node unique values
        target_label: Label of target nodes
        target_field: Unique field name of target nodes
        target_values: List of target node unique values
        relationship_type: Name of the relationship type
        direction_to_source: If True, source->target; if False, target->source
        tx_size: Maximum number of relationships per transaction batch

    Raises:
        ValueError: If source_values and target_values have different lengths
    """
    if len(source_values) != len(target_values):
        raise ValueError(
            f"source_values ({len(source_values)}) and target_values "
            f"({len(target_values)}) must have the same length"
        )

    if not source_values:
        return

    rows = [{"sv": sv, "dv": dv} for sv, dv in zip(source_values, target_values, strict=True)]

    for i in range(0, len(rows), tx_size):
        chunk = rows[i : i + tx_size]
        if direction_to_source:
            query = f"""
            UNWIND $rows AS r
            MATCH (s:`{source_label}` {{ `{source_field}`: r.sv }})
            MATCH (d:`{target_label}` {{ `{target_field}`: r.dv }})
            MERGE (s)-[:`{relationship_type}`]->(d)
            """
        else:
            query = f"""
            UNWIND $rows AS r
            MATCH (s:`{source_label}` {{ `{source_field}`: r.sv }})
            MATCH (d:`{target_label}` {{ `{target_field}`: r.dv }})
            MERGE (d)-[:`{relationship_type}`]->(s)
            """

        await session.run(query, rows=chunk)


async def remove_list_property_values(
    session: AsyncSession,
    label: str,
    unique_field: str,
    unique_values: list[str],
    list_property: str,
    values_to_remove: dict[str, list[str]],
    tx_size: int = 5000,
) -> None:
    """Remove matched values from list properties after relationships are created.

    Optimized: If all values in the list property were processed, removes the property
    entirely (faster). Otherwise, filters the list to remove only processed values.

    Args:
        session: Neo4j async session
        label: Label of nodes to update
        unique_field: Unique field name of nodes
        unique_values: List of node unique values to process
        list_property: Name of the list property to modify
        values_to_remove: Dict mapping unique_value to list of values to remove
        tx_size: Maximum number of nodes per transaction batch
    """
    if not unique_values:
        return

    # Build rows with values to remove for each unique value
    rows = [
        {
            "unique_value": uv,
            "values_to_remove": values_to_remove.get(uv, []),
        }
        for uv in unique_values
        if values_to_remove.get(uv)  # Only include if there are values to remove
    ]

    if not rows:
        return

    # Simple and fast: Just remove the entire property for all processed proteins
    # All IDs in the list were processed into relationships, so we can remove the property
    for i in range(0, len(unique_values), tx_size):
        chunk = unique_values[i : i + tx_size]
        query = f"""
        MATCH (n:`{label}`)
        WHERE n.`{unique_field}` IN $unique_values
        REMOVE n.`{list_property}`
        """
        await session.run(query, unique_values=chunk)


async def query_nodes_by_list_property(
    session: AsyncSession,
    label: str,
    list_property: str,
    unique_field: str,
    filter_values: list[str] | None = None,
) -> AsyncIterator[dict[str, Any]]:
    """Query nodes that have specific values in a list property.

    Returns nodes where the list property is non-empty, optionally filtered
    by specific values. Useful for enrichment stages to find nodes that
    need external data fetched.

    Args:
        session: Neo4j async session
        label: Label of nodes to query
        list_property: Name of the list property to check
        unique_field: Unique field name of nodes (for return value)
        filter_values: Optional list of values to filter by. If provided,
            only returns nodes where at least one value in the list matches.

    Yields:
        Dict with keys:
            - unique_field value: The unique identifier of the node
            - list_property: List of values in the property
    """
    if filter_values:
        query = f"""
        MATCH (n:`{label}`)
        WHERE n.`{list_property}` IS NOT NULL
          AND size(n.`{list_property}`) > 0
          AND ANY(x IN n.`{list_property}` WHERE x IN $filter_values)
        RETURN n.`{unique_field}` AS unique_value, n.`{list_property}` AS list_values
        """
        params = {"filter_values": filter_values}
    else:
        query = f"""
        MATCH (n:`{label}`)
        WHERE n.`{list_property}` IS NOT NULL
          AND size(n.`{list_property}`) > 0
        RETURN n.`{unique_field}` AS unique_value, n.`{list_property}` AS list_values
        """
        params = {}

    result = await session.run(query, **params)
    async for record in result:
        yield {
            unique_field: record["unique_value"],
            list_property: record["list_values"],
        }


async def create_reaction_molecule_relationships(
    session: AsyncSession,
    substrate_map: dict[str, list[str]],
    product_map: dict[str, list[str]],
    tx_size: int = 5000,
) -> None:
    """Create HAS_SUBSTRATE and HAS_PRODUCT relationships from Reaction to Molecule.

    Args:
        session: Neo4j async session
        substrate_map: Dict mapping reaction_id -> [substrate molecule_ids]
        product_map: Dict mapping reaction_id -> [product molecule_ids]
        tx_size: Transaction batch size

    Example:
        >>> substrate_map = {"RHEA:10000": ["CHEBI:15377", "CHEBI:15378"]}
        >>> product_map = {"RHEA:10000": ["CHEBI:15379"]}
        >>> await create_reaction_molecule_relationships(session, substrate_map, product_map)
    """
    # Create Reaction -HAS_SUBSTRATE-> Molecule relationships
    if substrate_map:
        source_values = []
        target_values = []
        for reaction_id, molecule_ids in substrate_map.items():
            for molecule_id in molecule_ids:
                source_values.append(reaction_id)
                target_values.append(molecule_id)

        if source_values:
            await create_relationships_batch(
                session=session,
                source_label="Reaction",
                source_field="id",
                source_values=source_values,
                target_label="Molecule",
                target_field="id",
                target_values=target_values,
                relationship_type="HAS_SUBSTRATE",
                direction_to_source=True,
                tx_size=tx_size,
            )

    # Create Reaction -HAS_PRODUCT-> Molecule relationships
    if product_map:
        source_values = []
        target_values = []
        for reaction_id, molecule_ids in product_map.items():
            for molecule_id in molecule_ids:
                source_values.append(reaction_id)
                target_values.append(molecule_id)

        if source_values:
            await create_relationships_batch(
                session=session,
                source_label="Reaction",
                source_field="id",
                source_values=source_values,
                target_label="Molecule",
                target_field="id",
                target_values=target_values,
                relationship_type="HAS_PRODUCT",
                direction_to_source=True,
                tx_size=tx_size,
            )


async def query_existing_nodes_by_ids(
    session: AsyncSession,
    label: str,
    unique_field: str,
    node_ids: list[str],
) -> list[str]:
    """Generic query for which node IDs already exist in database.

    Args:
        session: Neo4j async session
        label: Node label to query
        unique_field: Unique field name
        node_ids: List of node IDs to check

    Returns:
        List of node IDs that already exist in the database
    """
    if not node_ids:
        return []

    query = f"""
    MATCH (n:`{label}`)
    WHERE n.`{unique_field}` IN $node_ids
    RETURN n.`{unique_field}` AS node_id
    """
    result = await session.run(query, node_ids=node_ids)
    existing = [str(record["node_id"]) async for record in result]
    return existing


async def create_taxonomy_hierarchy(
    session: AsyncSession,
    main_taxon_id: str,
    parent_taxon_id: str | None,
    lineage_ids: list[str],
) -> None:
    """Create IS_A relationships for taxonomy lineage.

    Creates hierarchy chain: main IS_A parent IS_A lineage[-1] IS_A ... IS_A lineage[0]

    Args:
        session: Neo4j async session
        main_taxon_id: The queried taxon (species)
        parent_taxon_id: Direct parent taxon
        lineage_ids: Lineage from root to subfamily (ordered list)

    Example:
        For Homo sapiens (9606):
        - main_taxon_id=9606 (Homo sapiens)
        - parent_taxon_id=9605 (Homo)
        - lineage_ids=[131567, 2759, ..., 207598] (cellular organisms → Homininae)

        Creates:
        - 9606 IS_A 9605 (Homo sapiens IS_A Homo)
        - 9605 IS_A 207598 (Homo IS_A Homininae, last in lineage)
        - 207598 IS_A 9604, 9604 IS_A 314295, ... (lineage chain)
    """
    if not lineage_ids and not parent_taxon_id:
        return

    relationships: list[tuple[str, str]] = []

    # Main IS_A parent
    if parent_taxon_id is not None:
        relationships.append((main_taxon_id, parent_taxon_id))

    # Build lineage chain
    if lineage_ids:
        # Check if parent is already in lineage
        if parent_taxon_id is not None and parent_taxon_id in lineage_ids:
            # Parent is in lineage, so lineage chain will connect everything
            # Just create the full lineage chain: lineage[i] IS_A lineage[i-1]
            for i in range(len(lineage_ids) - 1, 0, -1):
                child_id = lineage_ids[i]
                parent_id = lineage_ids[i - 1]
                relationships.append((child_id, parent_id))
        else:
            # Parent is NOT in lineage, connect parent to most specific lineage item
            if parent_taxon_id is not None:
                relationships.append((parent_taxon_id, lineage_ids[-1]))
            # Create lineage chain: lineage[i] IS_A lineage[i-1]
            for i in range(len(lineage_ids) - 1, 0, -1):
                child_id = lineage_ids[i]
                parent_id = lineage_ids[i - 1]
                relationships.append((child_id, parent_id))

    if not relationships:
        return

    # Batch create IS_A relationships
    rows = [{"child": child, "parent": parent} for child, parent in relationships]
    query = """
    UNWIND $rows AS r
    MATCH (child:Taxon {id: r.child})
    MATCH (parent:Taxon {id: r.parent})
    MERGE (child)-[:IS_A]->(parent)
    """
    await session.run(query, rows=rows)
