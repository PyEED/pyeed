import re
import asyncio
import nest_asyncio
from typing import List, Union

from pyeed.core.dnaregion import DNARegion
from pyeed.core.proteinregiontype import ProteinRegionType
from pyeed.fetchers import AbstractFetcher, LOGGER
from pyeed.core.organism import Organism
from pyeed.core.proteininfo import ProteinInfo

from pyeed.fetchers.uniprotrequester import UniprotRequester
from pyeed.fetchers.interprorequester import InterProRequester


class UniprotFetcher(AbstractFetcher):

    def __init__(
        self, foreign_id: Union[int, List[int]], email: str = None, api_key: str = None
    ):
        super().__init__(foreign_id)

        self.api_key: str = api_key
        if email is None:
            self.email: str = self.get_substitute_email()
        self.taxonomy_dicts: List[dict] = None
        nest_asyncio.apply()

    def get(self):
        entries = UniprotRequester(self.foreign_id).make_request()
        return self._get_uniprot_result_dict(entries)

    def _get_uniprot_result_dict(self, records: List[dict]) -> dict:
        uniprot_dict = {}
        for record in records:
            uniprot_dict[record["accession"].upper()] = record
        return uniprot_dict

    def get_interpro_entries(self) -> dict:

        if not isinstance(self.foreign_id, list):
            self.foreign_id = [self.foreign_id]

        entries = asyncio.run(InterProRequester(self.foreign_id).make_request())

        return self._get_interpro_result_dict(entries)

    def _get_interpro_result_dict(self, interpro_records: List[dict]) -> dict:
        interpro_dict = {}
        for record in interpro_records:
            try:
                interpro_dict[
                    record["results"][0]["proteins"][0]["accession"].upper()
                ] = record
            except KeyError:
                LOGGER.info(f"Could not find the accession for {record}")
        return interpro_dict

    def fetch(self) -> List["ProteinInfo"]:
        uniprot_records_dict = self.get()
        interpro_records_dict = self.get_interpro_entries()

        protein_infos = []
        for key in uniprot_records_dict.keys():
            protein_infos.append(
                self.map(
                    uniprot_records_dict[key], interpro_records_dict[key], ProteinInfo
                )
            )

        return protein_infos

    def map(
        self,
        uniprot: dict,
        interpro: dict,
        cls: "ProteinInfo",
    ):

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
        protein_info = cls(
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

    def map_interpro(
        self, interpro: dict, protein_info: "ProteinInfo"
    ) -> "ProteinInfo":

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


if __name__ == "__main__":
    fetcher = UniprotFetcher(["P51587"])
    protein_info = fetcher.fetch()
    print(protein_info)
