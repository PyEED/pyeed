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
    Region,
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
        self.add_regions(record, protein)
        self.add_catalytic_activity(record, protein)
        self.add_go(record, protein)

    def add_sites(self, record: dict[str, Any], protein: Protein) -> None:
        data_dict: dict[str, list[int]] = defaultdict(list)

        for feature in record.get("features", []):
            if feature["type"] == "BINDING":
                for position in range(int(feature["begin"]), int(feature["end"]) + 1):
                    data_dict[feature["ligand"]["name"] + "$binding"].append(position)
            elif feature["type"] == "ACT_SITE":
                for position in range(int(feature["begin"]), int(feature["end"]) + 1):
                    data_dict[feature["category"] + "$site"].append(position)

        for entry, positions in data_dict.items():
            if entry.split("$")[1] == "binding":
                annotation = Annotation.BINDING_SITE.value
            elif entry.split("$")[1] == "site":
                annotation = Annotation.ACTIVE_SITE.value

            site = Site(
                name=entry.split("$")[0],
                annotation=annotation,
            )
            site.save()

            protein.site.connect(site, {"positions": positions})

    def add_regions(self, record: dict[str, Any], protein: Protein) -> None:
        data_list: list[tuple[str, tuple[int, int]]] = []

        for feature in record.get("features", []):
            if feature["type"] == "HELIX":
                data_list.append(
                    (
                        feature["category"] + "$helix",
                        (int(feature["begin"]), int(feature["end"])),
                    )
                )
            elif feature["type"] == "STRAND":
                data_list.append(
                    (
                        feature["category"] + "$strand",
                        (int(feature["begin"]), int(feature["end"])),
                    )
                )
            elif feature["type"] == "TURN":
                data_list.append(
                    (
                        feature["category"] + "$turn",
                        (int(feature["begin"]), int(feature["end"])),
                    )
                )
            elif feature["type"] == "SIGNAL":
                data_list.append(
                    (
                        feature["category"] + "$signal",
                        (int(feature["begin"]), int(feature["end"])),
                    )
                )
            elif feature["type"] == "PROPEP":
                data_list.append(
                    (
                        feature["category"] + "$propep",
                        (int(feature["begin"]), int(feature["end"])),
                    )
                )

        for name, positions in data_list:
            if name.split("$")[1] == "helix":
                annotation = Annotation.ALPHAHELIX.value
            elif name.split("$")[1] == "strand":
                annotation = Annotation.BETASTRAND.value
            elif name.split("$")[1] == "turn":
                annotation = Annotation.TURN.value
            elif name.split("$")[1] == "signal":
                annotation = Annotation.SIGNAL.value
            elif name.split("$")[1] == "propep":
                annotation = Annotation.PROPEP.value

            region = Region(
                name=name,
                annotation=annotation,
            )
            region.save()

            protein.region.connect(region, {"start": positions[0], "end": positions[1]})

    def add_catalytic_activity(self, record: dict[str, Any], protein: Protein) -> None:
        try:
            for reference in record["comments"]:
                if reference["type"] == "CATALYTIC_ACTIVITY":
                    catalytic_annotation = CatalyticActivity.get_or_save(
                        catalytic_id=int(reference["id"])
                        if reference.get("id")
                        else None,
                        name=reference["reaction"]["name"],
                    )
                    protein.catalytic_annotation.connect(catalytic_annotation)

        except Exception as e:
            logger.error(
                f"Error saving catalytic activity for {protein.accession_id}: {e}"
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
