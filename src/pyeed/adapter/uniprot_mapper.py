import json
from collections import defaultdict
from typing import Any

from httpx import Response
from loguru import logger

from pyeed.adapter.primary_db_adapter import PrimaryDBMapper
from pyeed.model import (
    Annotation,
    Reaction,
    GOAnnotation,
    Organism,
    Protein,
    Site,
)


class UniprotToPyeed(PrimaryDBMapper):
    def add_to_db(self, response: Response) -> None:
        # Organism information
        records = self.parse_response(response)

        for record in records:
            taxonomy_id = record["organism"]["taxonomy"]
            organism = Organism.get_or_save(
                taxonomy_id=taxonomy_id,
                name=record["organism"]["names"][0]["value"],
            )

            try:
                ec_number = record["protein"]["recommendedName"]["ecNumber"][0]["value"]
            except KeyError:
                ec_number = None

            # Protein name
            protein_data = record.get("protein", {})
            name = protein_data.get("recommendedName", {}).get("fullName", {}).get(
                "value"
            ) or protein_data.get("submittedName", [{}])[0].get("fullName", {}).get(
                "value"
            )

            try:
                protein = Protein.get_or_save(
                    accession_id=record["accession"],
                    sequence=record["sequence"]["sequence"],
                    mol_weight=float(record["sequence"]["mass"]),
                    ec_number=ec_number,
                    name=name,
                    seq_length=len(record["sequence"]["sequence"]),
                )
            except KeyError as e:
                logger.warning(
                    f"Error during mapping of {record['accession']} to graph model: {e}"
                )
                return

            protein.organism.connect(organism)
            self.add_reaction(record, protein)

        self.add_sites(record, protein)
        self.add_go(record, protein)

    def add_sites(self, record: dict[str, Any], protein: Protein) -> None:
        ligand_dict: dict[str, list[int]] = defaultdict(list)

        for feature in record.get("features", []):
            if feature["type"] == "BINDING":
                for position in range(int(feature["begin"]), int(feature["end"]) + 1):
                    ligand_dict[feature["ligand"]["name"]].append(position)

        for ligand, positions in ligand_dict.items():
            site = Site(
                name=ligand,
                annotation=Annotation.BINDING_SITE.value,
            )
            site.save()

            protein.site.connect(site, {"positions": positions})

    def add_reaction(self, record: dict[str, Any], protein: Protein) -> None:
        for reference in record.get("comments", []):  # Safe retrieval with .get()
            if reference.get("type") == "CATALYTIC_ACTIVITY":
                name = reference.get("reaction", {}).get("name", "")
                rhea_id = None  # Default value

                for db_ref in reference.get("reaction", {}).get("dbReferences", []):
                    if db_ref.get("id", "").startswith("RHEA:"):
                        rhea_id = db_ref["id"]
                        break  # Stop after finding the first match

                # Ensure we have both a reaction name and an RHEA ID
                if not name or not rhea_id:
                    logger.warning(f"Skipping {protein.accession_id}: Missing reaction name or RHEA ID")
                    continue  # Move to the next reference

                try:
                    left_side, right_side = name.split("=")
                    left_list = [s.strip() for s in left_side.split(" + ")]
                    right_list = [s.strip() for s in right_side.split(" + ")]

                    catalytic_annotation = Reaction.get_or_save(
                        rhea_id=rhea_id,
                        reactants=left_list,
                        products=right_list,
                    )
                    protein.reaction.connect(catalytic_annotation)

                except Exception as parse_error:
                    logger.error(f"Error processing reaction for {protein.accession_id}: {parse_error}")
                    continue  # Continue processing next accession_id


    def add_go(self, record: dict[str, Any], protein: Protein) -> None:
        for reference in record["dbReferences"]:
            if reference["type"] == "GO":
                go_annotation = GOAnnotation.get_or_save(
                    go_id=reference["id"],
                    term=reference["properties"]["term"],
                )

                protein.go_annotation.connect(go_annotation)

    def parse_response(self, response: Response) -> Any:
        """Parse UniProt JSON response."""
        if not response.content:
            return []

        try:
            return response.json()
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse UniProt response: {e}")
            return []
