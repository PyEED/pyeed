"""ChEBI ontology ingestion from OBO Graph JSON file.

Reads ChEBI data from a local JSON file (OBO Graph format) and upserts
Molecule nodes and their relationships to Neo4j.

Example:
    >>> from pyeed.ingest.chebi_from_file import ingest_chebi_from_file
    >>> from pyeed.db.neo4j import get_async_driver
    >>> import asyncio
    >>>
    >>> driver = get_async_driver()
    >>> asyncio.run(ingest_chebi_from_file(driver, "/path/to/chebi.json"))
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import httpx
from loguru import logger
from neo4j import AsyncDriver

from pyeed.db.query_utils import execute_write
from pyeed.ingest.model import Molecule

__all__ = [
    "create_molecule_relationships",
    "ingest_chebi_from_file",
    "normalize_relationship_label",
    "upsert_molecules_from_file",
]


# --- Label normalization ---


async def normalize_relationship_label(
    pred: str,
    client: httpx.AsyncClient | None = None,
    cache: dict[str, str] | None = None,
) -> str:
    """Extract and normalize relationship label from OLS API or return as-is.

    Args:
        pred: Input string - either a URL to purl.obolibrary.org/obo/ or a plain string.
        client: Optional reusable httpx client for batching.
        cache: Optional cache dict for memoization across calls.

    Returns:
        If input is a purl.obolibrary.org/obo/ URL: label from OLS API in UPPER_SNAKE_CASE.
        If input is not a URL: the input string in UPPER_SNAKE_CASE.

    Raises:
        ValueError: If input is a URL but not from purl.obolibrary.org/obo/.
        httpx.HTTPError: If OLS API request fails.
    """
    # Check cache first
    if cache is not None and pred in cache:
        return cache[pred]

    is_url = pred.startswith(("http://", "https://"))

    if not is_url:
        result = pred.upper().replace(" ", "_")
        if cache is not None:
            cache[pred] = result
        return result

    # Validate it's a purl.obolibrary.org/obo/ URL
    if "purl.obolibrary.org/obo/" not in pred:
        raise ValueError(f"Expected URL from purl.obolibrary.org/obo/, got: {pred}")

    # Fetch label from OLS API
    ols_url = "https://www.ebi.ac.uk/ols4/api/ontologies/ro/properties"
    params = {"iri": pred}

    if client is not None:
        resp = await client.get(ols_url, params=params)
    else:
        async with httpx.AsyncClient(timeout=30.0) as temp_client:
            resp = await temp_client.get(ols_url, params=params)

    resp.raise_for_status()
    data = resp.json()

    # Extract label from response
    properties = data.get("_embedded", {}).get("properties", [])
    if not properties:
        raise ValueError(f"No properties found in OLS response for {pred}")

    label = properties[0].get("label")
    if not label:
        raise ValueError(f"No label found in OLS response for {pred}")

    # Convert to UPPER_SNAKE_CASE
    result = label.upper().replace(" ", "_")
    if cache is not None:
        cache[pred] = result
    return result


def _extract_chebi_id(obo_iri: str) -> str:
    """Extract ChEBI ID from OBO IRI.

    Args:
        obo_iri: Full IRI like 'http://purl.obolibrary.org/obo/CHEBI_100862'

    Returns:
        ChEBI ID like 'CHEBI:100862'
    """
    local = obo_iri.rsplit("/", 1)[-1]  # "CHEBI_100862"
    return local.replace("_", ":")  # "CHEBI:100862"


# --- File reading ---


def read_chebi_json(path: str | Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Read ChEBI OBO Graph JSON and extract nodes and edges.

    Args:
        path: Path to ChEBI JSON file.

    Returns:
        Tuple of (nodes, edges) lists.

    Raises:
        FileNotFoundError: If file doesn't exist.
        json.JSONDecodeError: If file is not valid JSON.
        KeyError: If expected structure is missing.
    """
    path = Path(path)
    logger.info(f"Reading ChEBI ontology from {path}")

    with path.open() as f:
        data = json.load(f)

    graphs = data["graphs"]
    if not graphs:
        raise ValueError("No graphs found in ChEBI JSON")

    graph = graphs[0]
    nodes = graph.get("nodes", [])
    edges = graph.get("edges", [])

    logger.info(f"Found {len(nodes)} nodes and {len(edges)} edges")
    return nodes, edges


def parse_molecules(nodes: list[dict[str, Any]]) -> list[Molecule]:
    """Parse raw ChEBI nodes into Molecule objects.

    Args:
        nodes: List of node dicts from ChEBI JSON.

    Returns:
        List of Molecule objects.
    """
    molecules = []
    for node in nodes:
        # Skip nodes without proper ChEBI ID
        node_id = node.get("id", "")
        if "CHEBI_" not in node_id:
            continue
        try:
            mol = Molecule.from_chebi_node(node)
            molecules.append(mol)
        except Exception as e:
            logger.warning(f"Failed to parse node {node_id}: {e}")
    return molecules


# --- Relationship creation ---


async def create_molecule_relationships(
    driver: AsyncDriver,
    edges: list[dict[str, Any]],
    *,
    tx_size: int = 5000,
    resolve_labels: bool = True,
) -> dict[str, int]:
    """Create relationships between Molecule nodes from ChEBI edges.

    Args:
        driver: Neo4j async driver.
        edges: List of edge dicts from ChEBI JSON with 'sub', 'pred', 'obj'.
        tx_size: Transaction batch size.
        resolve_labels: If True, resolve OBO URLs to labels via OLS API.

    Returns:
        Dict mapping relationship type to count created.
    """
    # Group edges by predicate for batch processing
    edges_by_pred: dict[str, list[tuple[str, str]]] = defaultdict(list)

    for edge in edges:
        sub = edge.get("sub", "")
        pred = edge.get("pred", "")
        obj = edge.get("obj", "")

        # Only process ChEBI-to-ChEBI edges
        if "CHEBI_" not in sub or "CHEBI_" not in obj:
            continue

        sub_id = _extract_chebi_id(sub)
        obj_id = _extract_chebi_id(obj)
        edges_by_pred[pred].append((sub_id, obj_id))

    logger.info(f"Found {len(edges_by_pred)} unique relationship types")

    # Resolve relationship labels
    label_cache: dict[str, str] = {}
    rel_type_map: dict[str, str] = {}

    if resolve_labels:
        async with httpx.AsyncClient(timeout=30.0) as client:
            for pred in edges_by_pred:
                try:
                    rel_type = await normalize_relationship_label(pred, client, label_cache)
                    rel_type_map[pred] = rel_type
                    logger.debug(f"Resolved '{pred}' -> '{rel_type}'")
                except Exception as e:
                    logger.warning(f"Failed to resolve label for {pred}: {e}")
                    # Fallback: use last part of IRI or original string
                    fallback = pred.rsplit("/", 1)[-1] if "/" in pred else pred
                    rel_type_map[pred] = fallback.upper().replace(" ", "_")
    else:
        for pred in edges_by_pred:
            fallback = pred.rsplit("/", 1)[-1] if "/" in pred else pred
            rel_type_map[pred] = fallback.upper().replace(" ", "_")

    # Create relationships in batches
    counts: dict[str, int] = {}

    for pred, pairs in edges_by_pred.items():
        rel_type = rel_type_map[pred]
        logger.info(f"Creating {len(pairs)} '{rel_type}' relationships")

        cypher = f"""
        UNWIND $rows AS row
        MATCH (sub:Molecule {{id: row.sub_id}})
        MATCH (obj:Molecule {{id: row.obj_id}})
        MERGE (sub)-[r:`{rel_type}`]->(obj)
        """

        rows = [{"sub_id": sub_id, "obj_id": obj_id} for sub_id, obj_id in pairs]

        # Process in chunks
        for i in range(0, len(rows), tx_size):
            chunk = rows[i : i + tx_size]
            await execute_write(driver=driver, query=cypher, rows=chunk)

        counts[rel_type] = len(pairs)

    return counts


# --- Main ingestion functions ---


async def upsert_molecules_from_file(
    driver: AsyncDriver,
    path: str | Path,
    *,
    tx_size: int = 5000,
) -> int:
    """Read ChEBI JSON and upsert all Molecule nodes.

    Args:
        driver: Neo4j async driver.
        path: Path to ChEBI JSON file.
        tx_size: Transaction batch size.

    Returns:
        Number of molecules upserted.
    """
    nodes, _ = read_chebi_json(path)
    molecules = parse_molecules(nodes)

    logger.info(f"Upserting {len(molecules)} molecules")
    await Molecule._bulk_upsert(driver, molecules, tx_size=tx_size)

    return len(molecules)


async def ingest_chebi_from_file(
    driver: AsyncDriver,
    path: str | Path,
    *,
    tx_size: int = 5000,
    create_relationships: bool = True,
    resolve_labels: bool = True,
) -> dict[str, Any]:
    """Full ChEBI ingestion: upsert molecules and create relationships.

    Args:
        driver: Neo4j async driver.
        path: Path to ChEBI OBO Graph JSON file.
        tx_size: Transaction batch size for both nodes and relationships.
        create_relationships: If True, also create relationships between molecules.
        resolve_labels: If True, resolve OBO URLs to labels via OLS API.

    Returns:
        Dict with ingestion statistics.

    Example:
        >>> from pyeed.ingest.chebi_from_file import ingest_chebi_from_file
        >>> from pyeed.db.neo4j import get_async_driver
        >>> import asyncio
        >>>
        >>> driver = get_async_driver()
        >>> stats = asyncio.run(ingest_chebi_from_file(
        ...     driver,
        ...     "/path/to/chebi.json",
        ...     tx_size=5000,
        ... ))
        >>> print(stats)
        {'molecules': 224140, 'relationships': {'IS_A': 150000, ...}}
    """
    nodes, edges = read_chebi_json(path)

    # Step 1: Parse and upsert molecules
    molecules = parse_molecules(nodes)
    logger.info(f"Upserting {len(molecules)} molecules")
    await Molecule._bulk_upsert(driver, molecules, tx_size=tx_size)

    stats: dict[str, Any] = {"molecules": len(molecules)}

    # Step 2: Create relationships
    if create_relationships:
        rel_counts = await create_molecule_relationships(
            driver,
            edges,
            tx_size=tx_size,
            resolve_labels=resolve_labels,
        )
        stats["relationships"] = rel_counts
        stats["total_relationships"] = sum(rel_counts.values())

    logger.info(f"ChEBI ingestion complete: {stats}")
    return stats


# --- CLI ---

if __name__ == "__main__":
    import asyncio
    import sys

    from rich import print as rprint
    from rich.progress import Progress

    from pyeed.db.neo4j import get_async_driver

    async def main(path: str) -> None:
        driver = get_async_driver()

        with Progress() as progress:
            task = progress.add_task("[cyan]Ingesting ChEBI...", total=None)

            stats = await ingest_chebi_from_file(
                driver,
                path,
                tx_size=5000,
                create_relationships=True,
                resolve_labels=True,
            )

            progress.update(task, completed=True)

        rprint("[green]Done![/green]")
        rprint(stats)

    expected_args = 2  # script name + path
    if len(sys.argv) < expected_args:
        print("Usage: python -m pyeed.ingest.chebi_from_file <path_to_chebi.json>")
        sys.exit(1)

    asyncio.run(main(sys.argv[1]))
