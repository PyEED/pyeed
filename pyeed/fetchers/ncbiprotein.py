import re
from typing import List, Union
from pyeed.fetchers import AbstractFetcher, LOGGER, NCBITaxonomyFetcher
from pyeed.fetchers.entrezrequester import NCBIRequester

import Bio
from Bio.SeqRecord import SeqRecord
from pyeed.core.dnaregion import DNARegion
from pyeed.core.dnaregiontype import DNARegionType

from pyeed.core.organism import Organism
from pyeed.core.proteininfo import ProteinInfo
from pyeed.core.proteinregion import ProteinRegion
from pyeed.core.proteinsitetype import ProteinSiteType
from pyeed.core.span import Span


class NCBIProteinFetcher(AbstractFetcher):
    """
    The `NCBIProteinFetcher` class is a subclass of the `AbstractFetcher` class and is used
    to fetch protein data from the NCBI database.
    """

    def __init__(
        self, foreign_id: Union[int, List[int]], email: str = None, api_key: str = None
    ):
        super().__init__(foreign_id)

        self.api_key: str = api_key
        if email is None:
            self.email: str = self.get_substitute_email()
        self.taxonomy_dicts: List[dict] = None

    def get(self) -> List[SeqRecord]:
        """
        Fetches protein data from NCBI and returns a list of dictionaries of the results.
        """

        return NCBIRequester(
            foreign_id=self.foreign_id,
            email=self.email,
            db="protein",
            api_key=self.api_key,
            retmode="text",
            rettype="genbank",
        ).make_request()

    def map(
        self, seq_records: List[SeqRecord], cls: "ProteinInfo"
    ) -> List["ProteinInfo"]:
        """
        Maps the fetched data to an instance of the `ProteinInfo` class.
        """

        protein_infos = []
        for record in seq_records:
            protein_info = cls(source_id=record.id, sequence=str(record.seq))

            protein_info.organism = Organism(**self.map_organism(record))

            protein_info = self.map_protein(record, protein_info)

            protein_info = self.map_regions(record, protein_info)

            protein_info = self.map_sites(record, protein_info)

            protein_info = self.map_cds(record, protein_info)

            protein_infos.append(protein_info)

        return protein_infos

    def map_organism(self, seq_record: SeqRecord) -> dict:
        """
        Gets the organism name and taxonomy ID from the source data.
        Maps it to an Organism object.
        """

        feature = self.get_feature(seq_record, "source")
        if len(feature) != 1:
            LOGGER.debug(
                f"Multiple features ({len(feature)}) of type `source` found for {seq_record.id}: {feature}"
            )
        feature = feature[0]

        try:
            if len(feature.qualifiers["db_xref"]) != 1:
                LOGGER.info(
                    f"For {seq_record.id} {feature.qualifiers['db_xref']} taxonomy ID(s) were found, using the first one. Skipping organism assignment"
                )
                return None

            taxonomy_id = feature.qualifiers["db_xref"][0]

            if ":" in taxonomy_id:
                taxonomy_id = int(taxonomy_id.split(":")[1])

        except KeyError:
            LOGGER.debug(f"No taxonomy ID found for {seq_record.id}: {feature}")
            return None

        try:
            organism_name = feature.qualifiers["organism"]
        except KeyError:
            LOGGER.debug(
                f"No organism name found for {seq_record.id}: {feature[0].qualifiers}"
            )
            organism_name = None

        return {"name": organism_name[0], "taxonomy_id": taxonomy_id}

    def fetch(self, cls: "ProteinInfo") -> List[ProteinInfo]:
        """
        Fetches protein data from NCBI and returns a list of instances of the 'ProteinInfo' class.
        """

        seq_records = self.get()
        protein_infos = self.map(seq_records, cls)

        unique_tax_ids = set([info.organism.taxonomy_id for info in protein_infos])

        organisms = NCBITaxonomyFetcher(
            list(unique_tax_ids), self.email, self.api_key
        ).fetch(Organism)

        for protein_info in protein_infos:
            for organism in organisms:
                if protein_info.organism.taxonomy_id == organism.taxonomy_id:
                    protein_info.organism = organism

        return protein_infos

    def map_protein(self, seq_record: SeqRecord, protein_info: "ProteinInfo"):
        """Maps protein data from a `Bio.SeqRecord` to a `ProteinInfo` object."""

        protein = self.get_feature(seq_record, "Protein")
        if len(protein) == 0:
            LOGGER.debug(
                f"No protein feature found for {seq_record.id}: {seq_record.features}"
            )

            return protein_info

        if len(protein) > 1:
            LOGGER.debug(
                f"Multiple features ({len(protein)}) of type `Protein` found for {seq_record.id}"
            )

        protein = protein[0]
        try:
            protein_info.name = protein.qualifiers["product"][0]
        except KeyError:
            LOGGER.debug(
                f"No protein name found for {seq_record.id}: {protein.qualifiers}"
            )
            try:
                protein_info.name = protein.qualifiers["name"][0]
            except KeyError:
                LOGGER.debug(
                    f"No protein name found for {seq_record.id}: {protein.qualifiers}"
                )
                protein_info.name = None

        try:
            protein_info.mol_weight = protein.qualifiers["calculated_mol_wt"][0]
        except KeyError:
            LOGGER.debug(
                f"No molecular weight found for {seq_record.id}: {protein.qualifiers}"
            )
            protein_info.mol_weight = None

        try:
            protein_info.ec_number = protein.qualifiers["EC_number"][0]
        except KeyError:
            LOGGER.debug(
                f"No EC number found for {seq_record.id}: {protein.qualifiers}"
            )
            protein_info.ec_number = None

        return protein_info

    def map_regions(self, seq_record: SeqRecord, protein_info: "ProteinInfo"):
        """Maps region data from a `Bio.SeqRecord` to a `ProteinInfo` object."""

        regions = self.get_feature(seq_record, "region")
        for region in regions:
            try:
                if "db_xref" not in region.qualifiers:
                    db_xref = None
                else:
                    db_xref = region.qualifiers["db_xref"][0]

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
                        cross_reference=db_xref,
                    )
                )
            except KeyError:
                LOGGER.debug(
                    f"Incomplete region data found for {seq_record.id}: {region.qualifiers}, skipping region"
                )

        return protein_info

    def map_sites(self, seq_record: SeqRecord, protein_info: ProteinInfo):
        """Maps site data from a `Bio.SeqRecord` to a `ProteinInfo` object."""

        sites = self.get_feature(seq_record, "site")
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
                LOGGER.debug(
                    f"Incomplete site data found for {seq_record.id}: {site.qualifiers}, skipping site"
                )

        return protein_info

    def map_cds(self, seq_record: SeqRecord, protein_info: ProteinInfo):
        """Maps coding sequence data from a `Bio.SeqRecord` to a `ProteinInfo` object."""

        cds = self.get_feature(seq_record, "CDS")
        if len(cds) > 1:
            LOGGER.info(
                f"Multiple features ({len(cds)}) of type `CDS` found for {seq_record.id}"
            )

        try:
            cds = cds[0]
        except IndexError:
            LOGGER.debug(f"No CDS found for {seq_record.id}: {cds}")

            return protein_info

        try:
            protein_info.coding_sequence_ref = self.get_cds_regions(
                cds.qualifiers["coded_by"][0]
            )
        except IndexError:
            LOGGER.debug(
                f"No coding sequence reference found for {seq_record.id}: {cds.qualifiers}"
            )

        return protein_info

    @staticmethod
    def get_cds_regions(coded_by: dict) -> List[DNARegion]:
        """Extract coding sequence from the 'coded_by' qualifier and return a list of DNARegion objects."""

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

    def get_feature(
        self, seq_record: SeqRecord, feature_type: str
    ) -> "Bio.SeqFeature.SeqFeature":
        """Returns a list of features of a given type from a `Bio.SeqRecord` object."""
        return [
            feature
            for feature in seq_record.features
            if feature.type.lower() == feature_type.lower()
        ]


if __name__ == "__main__":
    id = ["UCS38941.1", "NP_001191"]
    fetcher = NCBIProteinFetcher(id)
    res = fetcher.fetch()
    print(res[0])
    print(res[1])
