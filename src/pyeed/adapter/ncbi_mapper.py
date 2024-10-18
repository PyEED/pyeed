import io
from abc import abstractmethod
from collections import defaultdict
from typing import Generic, TypeVar
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from loguru import logger

from pyeed.model import Annotation, GOAnnotation, Organism, Protein, Site

T = TypeVar("T")

class PrimaryDBtoPyeed(Generic[T]):
    @abstractmethod
    def add_to_db(self, data: dict):
        pass

class NCBIProteinToPyeed(PrimaryDBtoPyeed[Protein]):

    def get_feature(self, seq_record: SeqRecord, feature_type: str) -> SeqFeature:
        """Returns a list of features of a given type from a `Bio.SeqRecord` object."""
        return [
            feature
            for feature in seq_record.features
            if feature.type.lower() == feature_type.lower()
        ]

    def map_organism(self, seq_record: SeqRecord) -> dict:
        """
        Gets the organism name and taxonomy ID from the source data.
        Maps it to an Organism object.
        """

        feature = self.get_feature(seq_record, "source")
        if len(feature) != 1:
            logger.debug(
                f"Multiple features ({len(feature)}) of type `source` found for {seq_record.id}: {feature}"
            )
        feature = feature[0]

        try:
            if len(feature.qualifiers["db_xref"]) != 1:
                logger.info(
                    f"For {seq_record.id} {feature.qualifiers['db_xref']} taxonomy ID(s) were found, using the first one. Skipping organism assignment"
                )
                return None

            taxonomy_id = feature.qualifiers["db_xref"][0]

            if ":" in taxonomy_id:
                taxonomy_id = int(taxonomy_id.split(":")[1])

        except KeyError:
            logger.debug(f"No taxonomy ID found for {seq_record.id}: {feature}")
            return None

        try:
            organism_name = feature.qualifiers["organism"]
        except KeyError:
            logger.debug(
                f"No organism name found for {seq_record.id}: {feature[0].qualifiers}"
            )
            organism_name = None

        return organism_name[0], taxonomy_id

    def map_protein(self, seq_record: SeqRecord):

        protein_info_dict = {'name': None, 'mol_weight': None, 'ec_number': None}

        protein = self.get_feature(seq_record, "Protein")
        if len(protein) == 0:
            logger.debug(
                f"No protein feature found for {seq_record.id}: {seq_record.features}"
            )

            return protein_info_dict

        if len(protein) > 1:
            logger.debug(
                f"Multiple features ({len(protein)}) of type `Protein` found for {seq_record.id}"
            )

        protein = protein[0]
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
            protein_info_dict["mol_weight"] = protein.qualifiers["calculated_mol_wt"][0]
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

    """
    def map_regions(self, seq_record: SeqRecord, protein_info: ProteinRecord):

        regions = self.get_feature(seq_record, "region")
        for region in regions:
            try:
                if "db_xref" not in region.qualifiers:
                    db_xref = None
                else:
                    db_xref = region.qualifiers["db_xref"][0]

                protein_info.regions.append(
                    Region(
                        id=region.qualifiers["region_name"][0],
                        start=int(region.location.start),
                        end=int(region.location.end),
                        note=region.qualifiers["note"][0],
                        cross_reference=db_xref,
                    )
                )
            except KeyError:
                logger.debug(
                    f"Incomplete region data found for {seq_record.id}: {region.qualifiers}, skipping region"
                )

        return protein_info
    
    def map_sites(self, seq_record: SeqRecord, protein_info: ProteinRecord):

        sites = self.get_feature(seq_record, "site")
        for site in sites:
            try:
                protein_info.add_to_sites(
                    name=site.qualifiers["site_type"][0],
                    id=site.qualifiers["site_type"][0].lower(),
                    positions=[int(part.start) for part in site.location.parts],
                    cross_ref=site.qualifiers["db_xref"][0],
                )
            except KeyError:
                LOGGER.debug(
                    f"Incomplete site data found for {seq_record.id}: {site.qualifiers}, skipping site"
                )

        return protein_info

    def map_cds(self, seq_record: SeqRecord, protein_record: ProteinRecord):

        cds = self.get_feature(seq_record, "CDS")
        if len(cds) > 1:
            LOGGER.info(
                f"Multiple features ({len(cds)}) of type `CDS` found for {seq_record.id}"
            )

        try:
            cds = cds[0]
        except IndexError:
            LOGGER.debug(f"No CDS found for {seq_record.id}: {cds}")

            return protein_record

        try:
            protein_record.coding_sequence = self.get_cds_regions(
                cds.qualifiers["coded_by"][0]
            )
        except IndexError:
            LOGGER.debug(
                f"No coding sequence reference found for {seq_record.id}: {cds.qualifiers}"
            )

        return protein_record

    @staticmethod
    def get_cds_regions(coded_by: dict) -> List[DNARecord]:

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

        # TODO MAX ist hier nicht früher ein Fehler gewesen? es wird doch nur die etzte Region return, macht doch kein Sinn
        # Extract the start and end position of each region
        cds_ranges = [region.split(":")[1] for region in cds_regions]
        regions = []
        for region in cds_ranges:
            start, end = region.split("..")

            region = Region(
                id=reference_ids[0],
                start=int(start),
                end=int(end),
                type=Annotation.CODING_SEQ,
            )
            regions.append(region)

        return regions

    """

    def add_to_db(self, record: SeqIO.SeqRecord):

        logger.info(f"Mapping {record.id} to PyEED model")

        organism_name, taxonomy_id = self.map_organism(record)

        # Organism information
        organism = Organism.get_or_save(
            taxonomy_id=taxonomy_id,
            name=
            organism_name,
        )

        # Protein information
        protein_info_dict = self.map_protein(record)

        try:
            ec_number = protein_info_dict["ec_number"]
        except KeyError:
            ec_number = None

        try:
            protein = Protein.get_or_save(
                accession_id=record.id,
                sequence=record.seq.__str__(),
                mol_weight=protein_info_dict["mol_weight"],
                ec_number=ec_number,
                name=protein_info_dict["name"],
                seq_length=len(record.seq),
            )
        except KeyError as e:
            logger.warning(
                f"Error during mapping of {record.id} to graph model: {e}"
            )
            return

        protein.organism.connect(organism)

        # self.add_sites(data, protein)
        # self.add_go(data, protein)





    

