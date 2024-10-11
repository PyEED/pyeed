from abc import abstractmethod
from collections import defaultdict
from typing import Any, Generic, TypeVar

from pyeed.model import Annotation, Organism, ProteinRecord, Site

T = TypeVar("T")


class PrimaryDBResponseMapper(Generic[T]):
    @abstractmethod
    async def transform(self, data: dict) -> Any:
        pass


class UniprotMapper(PrimaryDBResponseMapper[ProteinRecord]):
    async def transform(self, data: dict) -> ProteinRecord:
        # Organism information
        organism = Organism(
            taxonomy_id=data["organism"]["taxonomy"],
        )
        print(organism)

        try:
            ec_number = data["protein"]["recommendedName"]["ecNumber"][0]["value"]
        except KeyError:
            ec_number = None

        sites = self.map_sites(data)

        protein_info = ProteinRecord(
            id=data["accession"],
            sequence=data["sequence"]["sequence"],
            name=data["protein"]["recommendedName"]["fullName"]["value"],
            ec_number=ec_number,
            mol_weight=data["sequence"]["mass"],
            organism=organism,
            sites=sites,
        )

        return protein_info

    def map_sites(self, data: dict) -> list[Site]:
        sites = []
        ligand_dict = defaultdict(list)
        for feature in data["features"]:
            if feature["type"] == "BINDING":
                if not feature["begin"] == feature["end"]:
                    continue
                ligand_dict[feature["ligand"]["name"]].append(int(feature["begin"]))

        for ligand, positions in ligand_dict.items():
            sites.append(
                Site(
                    name=ligand,
                    positions=positions,
                    annotation=Annotation.BINDING_SITE,
                )
            )

        return sites
