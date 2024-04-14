import re
import logging

from pyeed.core.dnaregion import DNARegion
from pyeed.core import ProteinInfo, Organism
from pyeed.core.proteinregiontype import ProteinRegionType


LOGGER = logging.getLogger(__name__)


class UniprotMapper:
    def __init__(self):
        pass

    def map(
        self,
        uniprot: dict,
        interpro: dict,
    ) -> ProteinInfo:
        """Maps the sequence information from Uniprot and annotations from InterPro
        records to a ProteinInfo object."""

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

        protein_info = ProteinInfo(
            source_id=uniprot["accession"],
            sequence=uniprot["sequence"]["sequence"],
            name=uniprot["protein"]["recommendedName"]["fullName"]["value"],
            ec_number=ec_number,
            mol_weight=uniprot["sequence"]["mass"],
            organism=organism,
        )
        for reference in uniprot["dbReferences"]:
            if reference["type"].upper() == "REFSEQ":
                try:
                    protein_info.coding_sequence_ref = DNARegion(
                        cross_reference=reference["properties"][
                            "nucleotide sequence ID"
                        ],
                    )
                except KeyError:
                    LOGGER.debug(
                        f"Could not find the coding sequence reference for {protein_info.source_id}"
                    )

        protein_info = self.map_interpro(interpro, protein_info)

        return protein_info

    def map_interpro(self, interpro: dict, protein_info: ProteinInfo) -> ProteinInfo:
        """Maps the InterPro records to a ProteinInfo object."""

        interpro_pattern = re.compile(r"IPR\d{6}")
        # pfam_pattern = re.compile(r"PF\d{5}")
        # panther_pattern = re.compile(r"PTHR\d{5}")

        for annotation in interpro["results"]:
            if interpro_pattern.search(annotation["metadata"]["accession"]):
                region = protein_info.add_to_regions(
                    name=annotation["metadata"]["name"],
                    type=ProteinRegionType.DOMAIN,
                    cross_reference=annotation["metadata"]["accession"],
                )
                region.add_to_spans(
                    start=annotation["proteins"][0]["entry_protein_locations"][0][
                        "fragments"
                    ][0]["start"],
                    end=annotation["proteins"][0]["entry_protein_locations"][0][
                        "fragments"
                    ][0]["end"],
                )

        return protein_info
