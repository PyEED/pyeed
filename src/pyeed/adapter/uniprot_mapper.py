from collections import defaultdict
from typing import Any

from loguru import logger
import traceback

from pyeed.adapter.primary_db_adapter import PrimaryDBtoPyeed
from pyeed.model import Annotation, GOAnnotation, Organism, Protein, Site, CatalyticActivity


class UniprotToPyeed(PrimaryDBtoPyeed):
    def add_to_db(self, data: Any) -> None:
        # Organism information
        taxonomy_id = data["organism"]["taxonomy"]
        organism = Organism.get_or_save(
            taxonomy_id=taxonomy_id,
            name=data["organism"]["names"][0]["value"],
        )

        try:
            ec_number = data["protein"]["recommendedName"]["ecNumber"][0]["value"]
        except KeyError:
            ec_number = None

        # Protein name
        protein_data = data.get("protein", {})
        name = protein_data.get("recommendedName", {}).get("fullName", {}).get(
            "value"
        ) or protein_data.get("submittedName", [{}])[0].get("fullName", {}).get("value")

        try:
            protein = Protein.get_or_save(
                accession_id=data["accession"],
                sequence=data["sequence"]["sequence"],
                mol_weight=float(data["sequence"]["mass"]),
                ec_number=ec_number,
                name=name,
                seq_length=len(data["sequence"]["sequence"]),
            )
        except KeyError as e:
            logger.warning(
                f"Error during mapping of {data['accession']} to graph model: {e}"
            )
            return

        protein.organism.connect(organism)

        self.add_sites(data, protein)
        self.add_catalytic_activity(data, protein)
        self.add_go(data, protein)

    def add_sites(self, data: dict, protein: Protein):
        ligand_dict: dict[str, list[int]] = defaultdict(list)

        for feature in data.get("features", []):
            if feature["type"] == "BINDING":
                for position in range(int(feature["begin"]), int(feature["end"]) + 1):
                    ligand_dict[feature["ligand"]["name"]].append(position)

        for ligand, positions in ligand_dict.items():
            site = Site(
                name=ligand,
                annotation=Annotation.BINDING_SITE.value,
            ).save()

            protein.site.connect(site, {"positions": positions})
    
    def add_catalytic_activity(self, data: dict, protein: Protein):
        try:
            for reference in data['comments']:
                if reference['type'] == "CATALYTIC_ACTIVITY":
                    catalytic_annotation = CatalyticActivity.get_or_save(
                        catalytic_id = int(reference["id"]) if reference.get("id") else None,
                        name=reference["reaction"]["name"],
                    )
                    protein.catalytic_annotation.connect(catalytic_annotation)
                    
        except Exception as e:
            logger.error(f"Error saving catalytic activity for {protein.accession_id}: {e}")
                
    def add_go(self, data: dict, protein: Protein):
        for reference in data["dbReferences"]:
            if reference["type"] == "GO":
                go_annotation = GOAnnotation.get_or_save(
                    go_id=reference["id"],
                    term=reference["properties"]["term"],
                )

                protein.go_annotation.connect(go_annotation)
