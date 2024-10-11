from abc import abstractmethod
from collections import defaultdict
from typing import Generic, TypeVar

from pyeed.model import Annotation, Organism, Protein, Site

T = TypeVar("T")


class PrimaryDBtoPyeed(Generic[T]):
    @abstractmethod
    async def add(self, data: dict):
        pass


class UniprotToPyeed(PrimaryDBtoPyeed[Protein]):
    async def add(self, data: dict):
        # Organism information
        taxonomy_id = data["organism"]["taxonomy"]

        organism = Organism(taxonomy_id=taxonomy_id).save()

        try:
            ec_number = data["protein"]["recommendedName"]["ecNumber"][0]["value"]
        except KeyError:
            ec_number = None

        protein = Protein(
            accession_id=data["accession"],
            sequence=data["sequence"]["sequence"],
        )
        protein.mol_weight = float(data["sequence"]["mass"])
        protein.ec_number = ec_number
        protein.name = data["protein"]["recommendedName"]["fullName"]["value"]
        protein.seq_length = len(protein.sequence)

        protein.save()

        protein.organism.connect(organism)
        organism.protein.connect(protein)

        self.add_sites(data, protein)

    def add_sites(self, data: dict, protein: Protein):
        ligand_dict = defaultdict(list)
        for feature in data["features"]:
            if feature["type"] == "BINDING":
                if not feature["begin"] == feature["end"]:
                    continue
                ligand_dict[feature["ligand"]["name"]].append(int(feature["begin"]))

        for ligand, positions in ligand_dict.items():
            site = Site(
                name=ligand,
                positions=positions,
                annotation=Annotation.BINDING_SITE.value,
            ).save()

            protein.sites.connect(site)
