import re
import logging
import secrets
import logging.config
from pathlib import Path
from abc import ABC, abstractmethod
from typing import Any, List

import Bio
from pyeed.core.dnaregion import DNARegion
from pyeed.core.dnaregiontype import DNARegionType

from pyeed.core.organism import Organism
from pyeed.core.proteininfo import ProteinInfo
from pyeed.core.proteinregion import ProteinRegion
from pyeed.core.proteinsitetype import ProteinSiteType
from pyeed.core.span import Span


path_config = Path(__file__).parent.parent.parent / "logging.conf"
logging.config.fileConfig(path_config)
LOGGER = logging.getLogger("pyeed")


class AbstractFetcher(ABC):
    def __init__(self, foreign_id: str):
        super().__init__()
        self.foreign_id = foreign_id

    @abstractmethod
    def get(self):
        pass

    @abstractmethod
    def map(self, handle: Any, cls):
        pass

    @staticmethod
    def get_substitute_email() -> str:
        return f"{secrets.token_hex(8)}@gmail.com"

    @staticmethod
    def make_chunks(input_list: list, chunk_size: int = 100) -> List[list]:
        """
        Splits a list into chunks of a given size.
        """
        if input_list is None:
            raise ValueError("input_list cannot be None.")

        if not isinstance(input_list, list):
            raise TypeError("input_list must be a list")

        return [
            input_list[i : i + chunk_size]
            for i in range(0, len(input_list), chunk_size)
        ]


class NCBIProteinParser(AbstractFetcher):

    def map(self, cls: "ProteinInfo"):

        protein_info = cls(source_id=self.source.id, sequence=str(self.source.seq))

        protein_info.organism = self.map_organism()
        protein_info = self.map_protein(protein_info)
        protein_info = self.map_regions(protein_info)
        protein_info = self.map_sites(protein_info)
        protein_info = self.map_cds(protein_info)

        return protein_info

    def map_organism(self) -> Organism:
        """
        Gets the organism name and taxonomy ID from the source data.
        Maps it to an Organism object.
        """

        feature = self.get_feature("source")
        if len(feature) != 1:
            LOGGER.debug(
                f"Multiple features ({len(feature)}) of type `source` found for {self.source.id}: {feature}"
            )
        feature = feature[0]

        try:
            if len(feature.qualifiers["db_xref"]) != 1:
                LOGGER.info(
                    f"For {self.source.id} {feature.qualifiers['db_xref']} taxonomy ID(s) were found, using the first one. Skipping organism assignment"
                )
                return None

            taxonomy_id = feature.qualifiers["db_xref"][0]

            if ":" in taxonomy_id:
                taxonomy_id = int(taxonomy_id.split(":")[1])

        except KeyError:
            LOGGER.debug(f"No taxonomy ID found for {self.source.id}: {feature}")
            return None

        try:
            organism_name = feature.qualifiers["organism"]
        except KeyError:
            LOGGER.debug(
                f"No organism name found for {self.source.id}: {feature[0].qualifiers}"
            )
            organism_name = None

        return Organism(name=organism_name[0], taxonomy_id=taxonomy_id)

    def map_protein(self, protein_info: ProteinInfo):

        protein = self.get_feature("Protein")
        if len(protein) == 0:
            LOGGER.debug(
                f"No protein feature found for {self.source.id}: {self.source.features}"
            )

            return protein_info

        if len(protein) > 1:
            LOGGER.debug(
                f"Multiple features ({len(protein)}) of type `Protein` found for {self.source.id}"
            )

        protein = protein[0]
        try:
            protein_info.name = protein.qualifiers["product"][0]
        except KeyError:
            LOGGER.debug(
                f"No protein name found for {self.source.id}: {protein.qualifiers}"
            )
            try:
                protein_info.name = protein.qualifiers["name"][0]
            except KeyError:
                LOGGER.debug(
                    f"No protein name found for {self.source.id}: {protein.qualifiers}"
                )
                protein_info.name = None

        try:
            protein_info.mol_weight = protein.qualifiers["calculated_mol_wt"][0]
        except KeyError:
            LOGGER.debug(
                f"No molecular weight found for {self.source.id}: {protein.qualifiers}"
            )
            protein_info.mol_weight = None

        try:
            protein_info.ec_number = protein.qualifiers["EC_number"][0]
        except KeyError:
            LOGGER.debug(
                f"No EC number found for {self.source.id}: {protein.qualifiers}"
            )
            protein_info.ec_number = None

        return protein_info

    def map_regions(self, protein_info: ProteinInfo):

        regions = self.get_feature("region")
        for region in regions:
            try:
                protein_info.regions.append(
                    ProteinRegion(
                        name=region.qualifiers["region_name"][0],
                        spans=[
                            Span(
                                start=int(region.location.start),
                                end=int(region.location.end),
                            )
                        ],
                        note=region.qualifiers["note"][0],
                        cross_reference=region.qualifiers["db_xref"][0],
                    )
                )
            except KeyError:
                LOGGER.debug(
                    f"Incomplete region data found for {self.source.id}: {region.qualifiers}, skipping region"
                )

        return protein_info

    def map_sites(self, protein_info: ProteinInfo):

        sites = self.get_feature("site")
        for site in sites:
            try:
                protein_info.add_to_sites(
                    name=site.qualifiers["site_type"][0],
                    type=ProteinSiteType.match_string(
                        site.qualifiers["site_type"][0].lower()
                    ),
                    positions=[int(part.start) for part in site.location.parts],
                    cross_ref=site.qualifiers["db_xref"][0],
                )
            except KeyError:
                LOGGER.warning(
                    f"Incomplete site data found for {self.source.id}: {site.qualifiers}, skipping site"
                )

        return protein_info

    def map_cds(self, protein_info: ProteinInfo):

        cds = self.get_feature("CDS")
        if len(cds) > 1:
            LOGGER.info(
                f"Multiple features ({len(cds)}) of type `CDS` found for {self.source.id}"
            )

        try:
            cds = cds[0]
        except IndexError:
            LOGGER.debug(f"No CDS found for {self.source.id}: {cds}")

            return protein_info

        try:
            protein_info.coding_sequence_ref = self.get_cds_regions(
                cds.qualifiers["coded_by"][0]
            )
        except IndexError:
            LOGGER.debug(
                f"No coding sequence reference found for {self.source.id}: {cds.qualifiers}"
            )

        return protein_info

    @staticmethod
    def get_cds_regions(coded_by: dict) -> List[DNARegion]:
        """Extract information about the coding sequence from the 'coded_by' qualifier."""

        cds_pattern = r"\w+\.\d+:\d+\.\.\d+\s?\d+"

        # Extract all regions from the 'coded_by' qualifier
        coded_by = coded_by.replace(">", "")
        coded_by = coded_by.replace("<", "")
        cds_regions = re.findall(cds_pattern, coded_by)
        cds_regions = [region.replace(" ", "") for region in cds_regions]

        # Extract the reference id from the first region
        reference_ids = [region.split(":")[0] for region in cds_regions]
        if not all(
            [reference_id == reference_ids[0] for reference_id in reference_ids]
        ):
            LOGGER.warning(
                "Nucleotide sequence references are not identical: {reference_ids}"
            )

        # Extract the start and end position of each region
        cds_ranges = [region.split(":")[1] for region in cds_regions]
        for region in cds_ranges:
            start, end = region.split("..")
            span = Span(start=int(start), end=int(end))  # noqa: F821

            region = DNARegion(
                id=reference_ids[0],
                spans=[span],
                type=DNARegionType.CODING_SEQUENCE,
            )

        return region

    def get_feature(self, feature_type: str) -> "Bio.SeqFeature.SeqFeature":
        return [
            feature
            for feature in self.source.features
            if feature.type.lower() == feature_type.lower()
        ]


class UniProtParser(AbstractFetcher):

    def map():
        pass


if __name__ == "__main__":
    from pyeed.ncbi.seq_io import get_ncbi_entry, get_ncbi_taxonomy
    from pyeed.core import Organism

    entry = get_ncbi_taxonomy("311400")

    parser = NCBITaxonomyParser(entry[0])
    print(parser.map(Organism))
