import re
from typing import Any, List, Tuple, TypeVar

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from loguru import logger

from pyeed.adapter.primary_db_adapter import PrimaryDBtoPyeed
from pyeed.model import Annotation, Organism, Protein, Region, Site

T = TypeVar("T")


class NCBIProteinToPyeed(PrimaryDBtoPyeed):

    def get_feature(self, seq_record: SeqRecord, feature_type: str) -> List[SeqFeature]:
        """Returns a list of features of a given type from a `Bio.SeqRecord` object."""
        return [
            feature
            for feature in seq_record.features
            if feature.type.lower() == feature_type.lower()
        ]

    def map_organism(self, seq_record: SeqRecord) -> Tuple[Any, Any]:
        """
        Gets the organism name and taxonomy ID from the source data.
        Maps it to an Organism object.
        """

        feature_list = self.get_feature(seq_record, "source")
        if len(feature_list) != 1:
            logger.debug(
                f"Multiple features ({len(feature_list)}) of type `source` found for {seq_record.id}: {feature_list}"
            )
        feature = feature_list[0]

        try:
            if len(feature.qualifiers["db_xref"]) != 1:
                logger.info(
                    f"For {seq_record.id} {feature.qualifiers['db_xref']} taxonomy ID(s) were found, using the first one. Skipping organism assignment"
                )
                return (None, None)

            taxonomy_id = feature.qualifiers["db_xref"][0]

            if ":" in taxonomy_id:
                taxonomy_id = int(taxonomy_id.split(":")[1])

        except KeyError:
            logger.debug(f"No taxonomy ID found for {seq_record.id}: {feature}")
            return (None, None)

        try:
            organism_name = feature.qualifiers["organism"]
        except KeyError:
            logger.debug(
                f"No organism name found for {seq_record.id}: {feature_list[0].qualifiers}"
            )
            organism_name = None

        return organism_name[0], taxonomy_id

    def map_protein(self, seq_record: SeqRecord) -> dict:
        """Maps protein information from a `Bio.SeqRecord` object to a dictionary.

        Args:
            seq_record (SeqRecord): A `Bio.SeqRecord` object.

        Returns:
            dict: A dictionary containing the protein information.
        """

        protein_info_dict: dict = {}
        protein_list = self.get_feature(seq_record, "Protein")
        if len(protein_list) == 0:
            logger.debug(
                f"No protein feature found for {seq_record.id}: {seq_record.features}"
            )

            return protein_info_dict

        if len(protein_list) > 1:
            logger.debug(
                f"Multiple features ({len(protein_list)}) of type `Protein` found for {seq_record.id}"
            )

        protein = protein_list[0]
        try:
            protein_info_dict["name"] = protein.qualifiers["product"][0]
        except KeyError:
            logger.debug(
                f"No protein name found for {seq_record.id}: {protein.qualifiers}"
            )
            try:
                protein_info_dict["name"] = protein.qualifiers["name"][0]
            except KeyError:
                logger.debug(
                    f"No protein name found for {seq_record.id}: {protein.qualifiers}"
                )
                protein_info_dict["name"] = None

        try:
            protein_info_dict["mol_weight"] = float(
                protein.qualifiers["calculated_mol_wt"][0]
            )
        except KeyError:
            logger.debug(
                f"No molecular weight found for {seq_record.id}: {protein.qualifiers}"
            )
            protein_info_dict["mol_weight"] = None

        try:
            protein_info_dict["ec_number"] = protein.qualifiers["EC_number"][0]
        except KeyError:
            logger.debug(
                f"No EC number found for {seq_record.id}: {protein.qualifiers}"
            )
            protein_info_dict["ec_number"] = None

        return protein_info_dict

    def map_sites(self, seq_record: SeqRecord):
        sites_list = []

        sites = self.get_feature(seq_record, "site")
        for site in sites:
            sites_dict = {}
            try:
                sites_dict["type"] = site.qualifiers["site_type"][0]
                sites_dict["positions"] = [
                    int(part.start) for part in site.location.parts
                ]
                sites_dict["cross_ref"] = site.qualifiers["db_xref"][0]

            except KeyError:
                logger.debug(
                    f"Incomplete site data found for {seq_record.id}: {site.qualifiers}, skipping site"
                )

            sites_list.append(sites_dict)

        return sites_list

    def add_sites(self, sites_list: List[dict], protein: Protein):
        for site_dict in sites_list:
            site = Site.get_or_save(
                name=site_dict["type"],
                annotation=Annotation.PROTEIN.value,
            )

            protein.site.connect(site, {"positions": site_dict["positions"]})

    def map_cds(self, seq_record: SeqRecord):
        cds = self.get_feature(seq_record, "CDS")


        if len(cds) > 1:
            logger.info(
                f"Multiple features ({len(cds)}) of type `CDS` found for {seq_record.id}"
            )

        try:
            cds_feature = cds[0]
        except IndexError:
            logger.debug(f"No CDS found for {seq_record.id}: {cds}")

            return None

        try:

            if "coded_by" not in cds_feature.qualifiers:
                logger.debug(
                    f"No coding sequence reference found for {seq_record.id}: {cds_feature.qualifiers}"
                )
                return None

            logger.info(f"CDS qualifiers: {cds_feature.qualifiers}")

            return self.get_cds_regions(cds_feature)
        except IndexError:
            logger.debug(
                f"No coding sequence reference found for {seq_record.id}: {cds_feature.qualifiers}"
            )

        return None

    @staticmethod
    def get_cds_regions(seq_feature: SeqFeature) -> List[dict]:
        coded_by =  seq_feature.qualifiers["coded_by"][0]
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
            logger.warning(
                "Nucleotide sequence references are not identical: {reference_ids}"
            )

        # Extract the start and end position of each region
        cds_ranges = [region.split(":")[1] for region in cds_regions]
        regions = []
        for i, region in enumerate(cds_ranges):
            start, end = region.split("..")

            region = {}
            region["start"] = int(start)
            region["end"] = int(end)
            region["type"] = Annotation.CODING_SEQ.value
            region["id"] = reference_ids[0]

            regions.append(region)

        return regions

    def map_regions(self, seq_record: SeqRecord):
        regions = self.get_feature(seq_record, "region")
        regions_list = []

        for region in regions:
            try:
                if "db_xref" not in region.qualifiers:
                    db_xref = None
                else:
                    db_xref = region.qualifiers["db_xref"][0]

                regions_list.append(
                    {
                        "id": region.qualifiers["region_name"][0],
                        "start": int(region.location.start),
                        "end": int(region.location.end),
                        "type": region.qualifiers["note"][0],
                        "cross_reference": db_xref,
                    }
                )
            except KeyError:
                logger.debug(
                    f"Incomplete region data found for {seq_record.id}: {region.qualifiers}, skipping region"
                )

        return regions_list

    def add_regions(self, regions_list: List[dict], protein: Protein):
        for region_dict in regions_list:
            region = Region.get_or_save(
                region_id=region_dict["id"],
                annotation=Annotation.PROTEIN.value,
            )

            protein.region.connect(
                region, {"start": region_dict["start"], "end": region_dict["end"]}
            )

    def add_to_db(self, record: SeqIO.SeqRecord):
        logger.info(f"Mapping {record.id} to PyEED model")

        organism_name, taxonomy_id = self.map_organism(record)

        # Organism information
        organism = Organism.get_or_save(
            taxonomy_id=taxonomy_id,
            name=organism_name,
        )

        # Protein information
        protein_info_dict = self.map_protein(record)

        logger.info(f"Adding protein {record.id} to the database")

        try:
            protein_info = {
                "sequence": str(record.seq),
                "mol_weight": protein_info_dict.get("mol_weight"),
                "ec_number": protein_info_dict.get("ec_number"),
                "name": protein_info_dict.get("name"),
                "seq_length": len(record.seq),
                "accession_id": str(record.id),
            }

            # Try to get existing protein
            protein = Protein.nodes.get_or_none(accession_id=record.id)
            if not protein:
                logger.debug(f"Protein {record.id} already exists")

            if protein is not None:
                logger.info(f"Updating existing protein {record.id}")
                protein.sequence = protein_info["sequence"]
                protein.mol_weight = protein_info["mol_weight"]
                protein.ec_number = protein_info["ec_number"]
                protein.name = protein_info["name"]
                protein.seq_length = protein_info["seq_length"]
                protein.save()

            else:
                logger.info(f"Creating new protein {record.id}")
                # Create new protein
                protein = Protein(**protein_info)
                protein.save()

        except KeyError as e:
            logger.warning(f"Error during mapping of {record.id} to graph model: {e}")
            return

        protein.organism.connect(organism)

        # Add the sites
        sites_list = self.map_sites(record)
        self.add_sites(sites_list, protein)

        # Add the coding sequence
        cds_regions = self.map_cds(record)
        if cds_regions is not None:
            for region in cds_regions:
                # add the id to protein nucleotide_id
                protein.nucleotide_id = region["id"]
                protein.nucleotide_start = region["start"]
                protein.nucleotide_end = region["end"]
                protein.save()

        # Add the regions
        regions_list = self.map_regions(record)
        self.add_regions(regions_list, protein)

        # Here we add the GO annotations
        # self.add_go(data, protein)
