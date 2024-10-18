from abc import abstractmethod
from collections import defaultdict
from typing import Generic, TypeVar

from loguru import logger

from pyeed.model import Annotation, GOAnnotation, Organism, Protein, Site

T = TypeVar("T")

class PrimaryDBtoPyeed(Generic[T]):
    @abstractmethod
    def add_to_db(self, data: dict):
        pass

class NCBIProteinToPyeed(PrimaryDBtoPyeed[Protein]):
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
        self.add_go(data, protein)
