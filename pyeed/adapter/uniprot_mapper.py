from abc import abstractmethod
from collections import defaultdict
from typing import Generic, TypeVar

from pyeed.model import Annotation, GOAnnotation, Organism, Protein, Site

T = TypeVar("T")


class PrimaryDBtoPyeed(Generic[T]):
    @abstractmethod
    def add_to_db(self, data: dict):
        pass


class UniprotToPyeed(PrimaryDBtoPyeed[Protein]):
    def add_to_db(self, data: dict):
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

        protein = Protein.get_or_save(
            accession_id=data["accession"],
            sequence=data["sequence"]["sequence"],
            mol_weight=float(data["sequence"]["mass"]),
            ec_number=ec_number,
            name=data["protein"]["recommendedName"]["fullName"]["value"],
            seq_length=len(data["sequence"]["sequence"]),
        )

        protein.organism.connect(organism)
        organism.protein.connect(protein)

        self.add_sites(data, protein)
        self.add_go(data, protein)

    def add_sites(self, data: dict, protein: Protein):
        ligand_dict = defaultdict(list)

        for feature in data.get("features", []):
            if feature["type"] == "BINDING":
                for position in range(int(feature["begin"]), int(feature["end"]) + 1):
                    ligand_dict[feature["ligand"]["name"]].append(position)

        for ligand, positions in ligand_dict.items():
            site = Site(
                name=ligand,
                positions=positions,
                annotation=Annotation.BINDING_SITE.value,
            ).save()

            protein.site.connect(site)

    def add_go(self, data: dict, protein: Protein):
        for reference in data["dbReferences"]:
            if reference["type"] == "GO":
                go_annotation = GOAnnotation.get_or_save(
                    go_id=reference["id"],
                    term=reference["properties"]["term"],
                )

                protein.go_annotation.connect(go_annotation)
