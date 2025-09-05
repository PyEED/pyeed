import asyncio
import json
import logging
import re
from dataclasses import dataclass
from logging import getLogger
from typing import Any, AsyncIterator

import httpx

from ..model import AnnotationType, Organism, Protein

logger = getLogger(__name__)

if not logger.handlers:
    fh = logging.FileHandler("protein_fetcher.log", mode="a")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(fh)

    # critical: set the logger level so DEBUG/INFO are emitted
    logger.setLevel(logging.DEBUG)

    # optional: avoid duplicates if root is configured elsewhere
    logger.propagate = False


@dataclass
class FetchConfig:
    """Simple config for HTTP fetching."""

    url: str
    headers: dict[str, str] | None = None
    timeout: int = 30
    rate_limit: int = 10
    max_concurrent: int = 5

    def __post_init__(self) -> None:
        if self.headers is None:
            self.headers = {"Accept": "application/json"}


async def fetch_uniprot_records(
    accessions: list[str], config: FetchConfig
) -> AsyncIterator[dict[str, Any]]:
    """
    Simple async fetcher - yields raw JSON records from UniProt.

    Args:
        accessions: List of UniProt accession IDs
        config: HTTP configuration

    Yields:
        Raw JSON records from UniProt API
    """
    if not accessions:
        return

    # Build OR query for multiple accessions
    query = " OR ".join(f"accession:{acc}" for acc in accessions)
    params = {
        "query": f"({query})",
        "format": "json",
        "size": "500",
        "fields": "id,sequence,accession,ft_domain,ft_site,ft_binding,organism_id,ft_act_site,rhea,ec,go",
    }

    logger.info(f"Fetching {len(accessions)} proteins from UniProt")

    async with httpx.AsyncClient(timeout=config.timeout) as client:
        try:
            response = await client.get(
                config.url, params=params, headers=config.headers
            )
            response.raise_for_status()

            # Parse JSON response
            data = response.json()
            records = data.get("results", []) if isinstance(data, dict) else data

            # save to json
            with open("uniprot_response.json", "w") as f:
                json.dump(data, f, indent=2)

            logger.info(f"Retrieved {len(records)} protein records")

            for record in records:
                yield record

        except httpx.HTTPStatusError as e:
            logger.error(f"HTTP error {e.response.status_code}: {e.response.text}")
            raise
        except httpx.RequestError as e:
            logger.error(f"Request error: {e}")
            raise
        except json.JSONDecodeError as e:
            logger.error(f"JSON decode error: {e}")
            raise


def map_uniprot_to_protein(record: dict[str, Any]) -> Protein:
    """
    Pure function: UniProt JSON record -> Protein object.
    Easy to test and reason about.

    Args:
        record: Raw UniProt JSON record

    Returns:
        Protein object
    """
    try:
        # Extract basic info
        accession = record["primaryAccession"]
        sequence_data = record["sequence"]
        sequence = sequence_data["value"]
        mol_weight = sequence_data["molWeight"]

        # Extract protein name
        protein_desc = record.get("proteinDescription", {})
        recommended_name = protein_desc.get("recommendedName", {})
        name = recommended_name.get("fullName", {}).get("value")

        # Try alternative names if no recommended name
        if not name:
            alt_names = protein_desc.get("alternativeNames", [])
            if alt_names:
                name = alt_names[0].get("fullName", {}).get("value")

        # Extract EC number
        ec_number = None
        ec_numbers = recommended_name.get("ecNumbers", [])
        if ec_numbers:
            ec_number = ec_numbers[0].get("value")

        # Create protein
        protein = Protein(
            accession_id=accession,
            sequence=sequence,
            mol_weight=mol_weight,
            ec_number=ec_number,
            name=name,
            seq_length=len(sequence),
            nucleotide_id=None,
            nucleotide_start=None,
            nucleotide_end=None,
            locus_tag=None,
        )

        # Add organism
        organism_data = record["organism"]
        organism = Organism(
            tax_id=organism_data["taxonId"],
            name=None,
        )
        protein.add_organism(protein.accession_id, organism)

        # Add features as annotations (robust)
        _add_features_robust(record, protein)

        # Add GO annotations (robust)
        _add_go_annotations(record, protein)

        # Add reactions (robust)
        _add_reactions(record, protein)

        return protein

    except KeyError as e:
        raise ValueError(f"Missing required field in UniProt record: {e}")
    except Exception as e:
        raise ValueError(f"Error mapping UniProt record: {e}")


def _add_reactions(record: dict[str, Any], protein: Protein) -> None:
    """
    Robustly extract UniProt reactions.
    Handles missing fields gracefully without raising exceptions.
    Only adds Rhea reaction IDs (not substrates or products).
    """
    from ..model import Reaction

    comments = record.get("comments", [])
    rhea_pattern = re.compile(r"^RHEA:\d+$")
    added_count = 0
    skipped_count = 0

    for comment in comments:
        if comment.get("commentType", "").lower().strip() != "catalytic activity":
            continue
        reaction = comment.get("reaction", {})
        if not reaction:
            logger.debug(f"Skipping invalid reaction: {reaction}")
            skipped_count += 1
            continue

        # Extract Rhea IDs from reactionCrossReferences
        cross_refs = reaction.get("reactionCrossReferences", [])
        for ref in cross_refs:
            db = ref.get("database", "")
            rid = ref.get("id", "")
            if db == "Rhea" and rhea_pattern.match(rid):
                # Only add the Rhea ID as a reaction, not substrates or products
                reaction_obj = Reaction(
                    rhea_id=rid,
                    description=reaction.get("name", None),
                )
                protein.add_reaction(protein.accession_id, reaction_obj)
                added_count += 1

    if added_count == 0:
        logger.debug(f"No Rhea reactions found for protein {protein.accession_id}")


def _add_go_annotations(record: dict[str, Any], protein: Protein) -> None:
    """
    Robustly extract UniProt GO annotations.
    Handles missing fields gracefully without raising exceptions.
    """
    from ..model import GOAnnotation

    cross_references = record.get("uniProtKBCrossReferences", [])
    if not cross_references:
        logger.debug(f"No cross references found for protein {protein.accession_id}")
        return

    processed_count = 0
    skipped_count = 0

    for cross_ref in cross_references:
        try:
            # Only process GO terms
            database = cross_ref.get("database", "")
            if database != "GO":
                continue

            # Extract GO ID
            go_id = cross_ref.get("id", "")
            if not go_id or not go_id.startswith("GO:"):
                logger.debug(f"Skipping invalid GO ID: {go_id}")
                skipped_count += 1
                continue

            # Extract GO term and category from properties
            properties = cross_ref.get("properties", [])
            go_term = None
            evidence_type = None

            for prop in properties:
                key = prop.get("key", "")
                value = prop.get("value", "")

                if key == "GoTerm":
                    go_term = value
                elif key == "GoEvidenceType":
                    evidence_type = value

            if not go_term:
                logger.debug(f"Skipping GO term {go_id}: no GoTerm property")
                skipped_count += 1
                continue

            # Create GO annotation
            go_annotation = GOAnnotation(
                go_id=go_id,
                term=go_term,
                definition=None,
            )

            # Add to protein
            protein.add_go_annotation(protein.accession_id, go_annotation)
            processed_count += 1

        except Exception as e:
            logger.warning(
                f"Error processing GO annotation for {protein.accession_id}: {e}"
            )
            skipped_count += 1
            continue

    logger.debug(
        f"GO annotation processing for {protein.accession_id}: "
        f"{processed_count} added, {skipped_count} skipped"
    )


def _add_features_robust(record: dict[str, Any], protein: Protein) -> None:
    """
    Robustly extract UniProt features as annotations.
    Handles missing fields gracefully without raising exceptions.
    """
    from ..model import Annotation

    features = record.get("features", [])
    if not features:
        logger.debug(f"No features found for protein {protein.accession_id}")
        return

    # Mapping from UniProt feature types to our annotation types
    feature_type_mapping = {
        "active site": AnnotationType.ACTIVE_SITE,
        "site": AnnotationType.SITE,
        "domain": AnnotationType.DOMAIN,
        "binding site": AnnotationType.BINDING_SITE,
        "metal ion-binding site": AnnotationType.BINDING_SITE,
        "nucleotide phosphate-binding region": AnnotationType.BINDING_SITE,
    }

    processed_count = 0
    skipped_count = 0

    for feature in features:
        try:
            # Extract feature type safely
            feature_type = feature.get("type", "").lower().strip()
            if not feature_type:
                skipped_count += 1
                continue

            # Check if we handle this feature type
            annotation_type = feature_type_mapping.get(feature_type)
            if not annotation_type:
                logger.debug(f"Skipping unsupported feature type: {feature_type}")
                skipped_count += 1
                continue

            # Extract location safely
            location = feature.get("location", {})
            if not location:
                logger.debug(f"Skipping feature {feature_type}: no location data")
                skipped_count += 1
                continue

            # Extract start and end positions
            start_data = location.get("start", {})
            end_data = location.get("end", {})

            start_val = start_data.get("value") if start_data else None
            end_val = end_data.get("value") if end_data else None

            if start_val is None or end_val is None:
                logger.debug(
                    f"Skipping feature {feature_type}: missing start/end values"
                )
                skipped_count += 1
                continue

            # Convert to integers safely
            try:
                start_pos = int(start_val)
                end_pos = int(end_val)
            except (ValueError, TypeError):
                logger.debug(
                    f"Skipping feature {feature_type}: invalid position values"
                )
                skipped_count += 1
                continue

            # Validate position range
            if start_pos <= 0 or end_pos <= 0 or start_pos > end_pos:
                logger.debug(f"Skipping feature {feature_type}: invalid position range")
                skipped_count += 1
                continue

            # Create position list
            positions = list(range(start_pos, end_pos + 1))

            # Extract additional metadata safely
            description = feature.get("description", None)
            if not description:
                description = None

            # Create annotation
            annotation = Annotation(
                annotation_type=annotation_type,
                positions=positions,
                description=description,
            )

            # Add to protein (accumulate all annotations)
            protein.add_sequence_annotation(protein.accession_id, annotation)
            processed_count += 1

        except Exception as e:
            logger.warning(f"Error processing feature for {protein.accession_id}: {e}")
            skipped_count += 1
            continue

    logger.debug(
        f"Feature processing for {protein.accession_id}: "
        f"{processed_count} added, {skipped_count} skipped"
    )


# Simple high-level interface
async def get_proteins_from_uniprot(accessions: list[str]) -> list[Protein]:
    """
    Get proteins from UniProt - simple and clean interface.

    Args:
        accessions: List of UniProt accession IDs

    Returns:
        List of Protein objects
    """
    config = FetchConfig(url="https://rest.uniprot.org/uniprotkb/search")

    proteins = []
    async for record in fetch_uniprot_records(accessions, config):
        try:
            protein = map_uniprot_to_protein(record)
            proteins.append(protein)
            logger.debug(f"Mapped protein {protein.accession_id}")
        except Exception as e:
            accession = record.get("primaryAccession", "unknown")
            logger.warning(f"Failed to map protein {accession}: {e}")

    logger.info(f"Successfully created {len(proteins)} protein objects")
    return proteins


# Optional buffer for batch operations
class ProteinBuffer:
    """Simple buffer for batch operations."""

    def __init__(self) -> None:
        self.proteins: list[Protein] = []

    def add(self, protein: Protein) -> None:
        """Add protein to buffer."""
        self.proteins.append(protein)

    def add_all(self, proteins: list[Protein]) -> None:
        """Add multiple proteins to buffer."""
        self.proteins.extend(proteins)

    def clear(self) -> None:
        """Clear buffer."""
        self.proteins.clear()

    def size(self) -> int:
        """Get buffer size."""
        return len(self.proteins)

    def flush_to_neo4j(self, writer: Any) -> int:
        """
        Batch write to Neo4j.

        Args:
            writer: Neo4j writer instance

        Returns:
            Number of proteins written
        """
        count = 0
        for protein in self.proteins:
            try:
                writer.create(protein)
                count += 1
            except Exception as e:
                logger.error(f"Failed to write protein {protein.accession_id}: {e}")

        self.clear()
        logger.info(f"Flushed {count} proteins to Neo4j")
        return count

    def save_to_file(self, path: str) -> None:
        """Save proteins as JSON file."""
        data = [p.model_dump() for p in self.proteins]
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
        logger.info(f"Saved {len(data)} proteins to {path}")


if __name__ == "__main__":
    from rich import print

    ids = ["P21802", "P69905", "P12345"]
    proteins = asyncio.run(get_proteins_from_uniprot(ids))

    print(proteins[-1])
