from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING

from pyeed.core import Organism
from pyeed.core.ontology import Ontology
from pyeed.core.region import Region

if TYPE_CHECKING:
    from pyeed.core import ProteinRecord

LOGGER = logging.getLogger(__name__)


class UniprotMapper:
    def __init__(self):
        pass

    def map_uniprot_data(self, uniprot_data: dict) -> ProteinRecord:
        from pyeed.core import ProteinRecord

        # Organism information
        organism = Organism(
            taxonomy_id=uniprot_data["organism"]["taxonomy"],
        )

        try:
            ec_number = uniprot_data["protein"]["recommendedName"]["ecNumber"][0][
                "value"
            ]
        except KeyError:
            ec_number = None

        protein_info = ProteinRecord(
            id=uniprot_data["accession"],
            sequence=uniprot_data["sequence"]["sequence"],
            name=uniprot_data["protein"]["recommendedName"]["fullName"]["value"],
            ec_number=ec_number,
            mol_weight=uniprot_data["sequence"]["mass"],
            organism=organism,
        )

        global_go_annotations = [
            Ontology.GO.value + go_annotation["id"]
            for go_annotation in uniprot_data["dbReferences"]
            if go_annotation["type"] == "GO"
        ]

        return protein_info

    def map(
        self,
        uniprot: dict,
        interpro: dict,
    ) -> ProteinRecord:
        """Maps the sequence information from Uniprot and annotations from InterPro
        records to a ProteinInfo object."""

        from pyeed.core import ProteinRecord

        assert (
            interpro["results"][0]["proteins"][0]["accession"].upper()
            == uniprot["accession"].upper()
        )

        organism = Organism(
            taxonomy_id=uniprot["organism"]["taxonomy"],
        )

        try:
            ec_number = uniprot["protein"]["recommendedName"]["ecNumber"][0]["value"]
        except KeyError:
            ec_number = None

        protein_info = ProteinRecord(
            id=uniprot["accession"],
            accession_id=uniprot["accession"],
            sequence=uniprot["sequence"]["sequence"],
            name=uniprot["protein"]["recommendedName"]["fullName"]["value"],
            ec_number=ec_number,
            mol_weight=uniprot["sequence"]["mass"],
            organism=organism,
        )
        for reference in uniprot["dbReferences"]:
            if reference["type"].upper() == "REFSEQ":
                try:
                    protein_info.coding_sequence.append(
                        Region(
                            id=reference["properties"]["nucleotide sequence ID"],
                        )
                    )
                except KeyError:
                    LOGGER.debug(
                        f"Could not find the coding sequence reference for {protein_info.id}"
                    )

        protein_info = self.map_interpro(interpro, protein_info)

        return protein_info

    def map_interpro(
        self, interpro: dict, protein_info: ProteinRecord
    ) -> ProteinRecord:
        """Maps the InterPro records to a ProteinInfo object."""

        interpro_pattern = re.compile(r"IPR\d{6}")
        # pfam_pattern = re.compile(r"PF\d{5}")
        # panther_pattern = re.compile(r"PTHR\d{5}")

        for annotation in interpro["results"]:
            if interpro_pattern.search(annotation["metadata"]["accession"]):
                region = protein_info.add_to_regions(
                    name=annotation["metadata"]["name"],
                    cross_reference=annotation["metadata"]["accession"],
                    start=annotation["proteins"][0]["entry_protein_locations"][0][
                        "fragments"
                    ][0]["start"],
                    end=annotation["proteins"][0]["entry_protein_locations"][0][
                        "fragments"
                    ][0]["end"],
                )

        return protein_info
