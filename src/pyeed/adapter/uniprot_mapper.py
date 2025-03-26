import json
from collections import defaultdict
from typing import Any

from httpx import Response
from loguru import logger

from pyeed.adapter.primary_db_adapter import PrimaryDBMapper
from pyeed.model import (
    Annotation,
    CatalyticActivity,
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

        self.add_sites(record, protein)
        self.add_catalytic_activity(record, protein)
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

    def add_catalytic_activity(self, record: dict[str, Any], protein: Protein) -> None:
        try:
            for reference in record["comments"]:
                if reference["type"] == "CATALYTIC_ACTIVITY":
                    name = reference["reaction"]["name"]
                    for i in reference["reaction"]["dbReferences"]:
                        if i['id'].startswith("RHEA:"):
                            rhea_id = i['id']
                            break
                    left_side, right_side = name.split("=")

                    # Further split each side by "+"
                    left_list = list(left_side.strip().split(" + "))
                    right_list = list(right_side.strip().split(" + "))
                    
                    catalytic_annotation = CatalyticActivity.get_or_save(
                        catalytic_id=int(reference["id"])
                        if reference.get("id")
                        else None,
                        rhea_id=rhea_id,
                        reactants = left_list,
                        products = right_list,
                    )
                    protein.catalytic_annotation.connect(catalytic_annotation)

        except Exception as e:
            logger.error(
                f"No Catalytic Activity for {protein.accession_id}: {e}"
            )

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
