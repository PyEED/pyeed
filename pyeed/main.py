import asyncio
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from neo4j import GraphDatabase

from .fetch.uniprot import get_proteins_from_uniprot_batched
from .model import Protein
from .neo4j_integration import PyeedNeo4jWriter

logger = logging.getLogger(__name__)


async def from_uniprot(
    ids: List[str],
    neo4j_uri: Optional[str] = None,
    neo4j_user: Optional[str] = None,
    neo4j_password: Optional[str] = None,
    output_file: Optional[str] = None,
    batch_size: int = 100,
    max_concurrent: int = 5,
    check_existing: bool = True,
) -> Dict[str, Any]:
    """
    High-performance method to fetch proteins from UniProt and optionally store in Neo4j.

    Args:
        ids: List of UniProt accession IDs
        neo4j_uri: Neo4j database URI (if None, saves to file)
        neo4j_user: Neo4j username
        neo4j_password: Neo4j password
        output_file: Output file path (default: proteins_{count}_records.json)
        batch_size: Number of IDs to process per batch
        max_concurrent: Maximum concurrent requests
        check_existing: Whether to check for existing proteins in Neo4j

    Returns:
        Dictionary with processing results
    """
    if all([neo4j_uri is None, neo4j_user is None, neo4j_password is None]):
        to_db = False
    else:
        to_db = True
    logger.info(f"Processing {len(ids)} UniProt IDs")

    # Determine IDs to fetch
    ids_to_fetch = set(ids)
    existing_count = 0

    # Check existing proteins in Neo4j if database provided
    if neo4j_uri and check_existing:
        logger.info("Checking for existing proteins in Neo4j...")
        existing_ids = await _check_existing_proteins(
            ids, neo4j_uri, neo4j_user, neo4j_password
        )
        existing_count = len(existing_ids)
        ids_to_fetch = ids_to_fetch - existing_ids
        logger.info(
            f"Found {existing_count} existing proteins, fetching {len(ids_to_fetch)} new ones"
        )

    # Fetch proteins from UniProt in batches
    proteins = []
    if ids_to_fetch:
        logger.info(f"Fetching {len(ids_to_fetch)} proteins from UniProt...")
        proteins = await get_proteins_from_uniprot_batched(
            list(ids_to_fetch), batch_size=batch_size, max_concurrent=max_concurrent
        )
        logger.info(f"Successfully fetched {len(proteins)} proteins")

    # Store results
    result = {
        "total_requested": len(ids),
        "existing_in_db": existing_count,
        "fetched_from_uniprot": len(proteins),
        "proteins_processed": len(proteins),
        "failed": len(ids_to_fetch) - len(proteins),
    }

    if proteins:
        if neo4j_uri:
            # Store in Neo4j
            logger.info(f"Storing {len(proteins)} proteins in Neo4j...")
            neo4j_results = await _store_proteins_neo4j(
                proteins, neo4j_uri, neo4j_user, neo4j_password, batch_size
            )
            result.update(neo4j_results)
        else:
            # Store in file
            output_path = output_file or f"proteins_{len(proteins)}_records.json"
            await _store_proteins_file(proteins, output_path)
            result["output_file"] = output_path
            logger.info(f"Saved {len(proteins)} proteins to {output_path}")

    return result


async def _check_existing_proteins(
    ids: List[str], uri: str, user: str, password: str
) -> Set[str]:
    """Check which protein IDs already exist in Neo4j."""
    existing_ids = set()

    try:
        driver = GraphDatabase.driver(uri, auth=(user, password))

        # Query in batches to avoid large parameter lists
        batch_size = 1000
        for i in range(0, len(ids), batch_size):
            batch_ids = ids[i : i + batch_size]

            with driver.session() as session:
                cypher = """
                UNWIND $ids as protein_id
                MATCH (p:Protein {accession_id: protein_id})
                RETURN p.accession_id as accession_id
                """
                result = session.run(cypher, {"ids": batch_ids})

                for record in result:
                    existing_ids.add(record["accession_id"])

        driver.close()

    except Exception as e:
        logger.warning(f"Could not check existing proteins: {e}")
        # Continue without checking - better to have duplicates than miss data

    return existing_ids


async def _store_proteins_neo4j(
    proteins: List[Protein], uri: str, user: str, password: str, batch_size: int = 100
) -> Dict[str, Any]:
    """Store proteins in Neo4j using batch transactions."""
    logger.info(
        f"Storing {len(proteins)} proteins in Neo4j with batch size {batch_size}"
    )

    writer = PyeedNeo4jWriter(uri, user, password)

    try:
        # Setup database constraints first
        writer.setup_database()

        # Process proteins in batches
        results = {"neo4j_success": 0, "neo4j_failed": 0, "neo4j_errors": []}

        for i in range(0, len(proteins), batch_size):
            batch = proteins[i : i + batch_size]
            logger.info(
                f"Processing Neo4j batch {i//batch_size + 1}/{(len(proteins)-1)//batch_size + 1}"
            )

            for protein in batch:
                try:
                    result = writer.add_protein(protein)
                    if result["success"]:
                        results["neo4j_success"] += 1
                    else:
                        results["neo4j_failed"] += 1
                        results["neo4j_errors"].append(
                            {
                                "protein_id": protein.accession_id,
                                "error": result.get("error", "Unknown error"),
                            }
                        )
                except Exception as e:
                    results["neo4j_failed"] += 1
                    results["neo4j_errors"].append(
                        {"protein_id": protein.accession_id, "error": str(e)}
                    )
                    logger.error(f"Failed to store protein {protein.accession_id}: {e}")

        logger.info(
            f"Neo4j storage complete: {results['neo4j_success']} success, {results['neo4j_failed']} failed"
        )
        return results

    finally:
        writer.close()


async def _store_proteins_file(proteins: List[Protein], output_path: str) -> None:
    """Store proteins as JSON file."""
    data = [protein.model_dump() for protein in proteins]

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        json.dump(
            {
                "metadata": {
                    "count": len(proteins),
                    "source": "UniProt",
                    "format": "Pyeed Protein objects",
                },
                "proteins": data,
            },
            f,
            indent=2,
        )


# Example usage and CLI interface
if __name__ == "__main__":
    import sys

    # Simple CLI interface
    if len(sys.argv) < 2:
        print(
            "Usage: python -m pyeed.main <ids...> [--neo4j-uri URI] [--output file.json]"
        )
        print(
            "Example: python -m pyeed.main P69905 P21802 --neo4j-uri bolt://localhost:7687"
        )
        print("Example: python -m pyeed.main P69905 P21802 --output my_proteins.json")
        sys.exit(1)

    # Parse arguments
    ids = []
    neo4j_uri = None
    output_file = None

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == "--neo4j-uri" and i + 1 < len(sys.argv):
            neo4j_uri = sys.argv[i + 1]
            i += 2
        elif arg == "--output" and i + 1 < len(sys.argv):
            output_file = sys.argv[i + 1]
            i += 2
        elif not arg.startswith("--"):
            ids.append(arg)
            i += 1
        else:
            i += 1

    # Run the main function
    async def main():
        result = await from_uniprot(
            ids=ids,
            neo4j_uri=neo4j_uri,
            neo4j_user="neo4j" if neo4j_uri else None,
            neo4j_password="12345678" if neo4j_uri else None,
            output_file=output_file,
        )
        print(f"Processing complete: {result}")

    asyncio.run(main())
