from typing import List, TypeVar

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from loguru import logger

from pyeed.adapter.primary_db_adapter import PrimaryDBtoPyeed
from pyeed.model import DNA, Annotation, Organism, Protein, Region, Site

T = TypeVar("T")


class NCBIDNAToPyeed(PrimaryDBtoPyeed):
    """
    Maps DNA sequence entries from NCBI to the PyEED graph object model and saves
    them to the database.
    """

    def add_to_db(self, record: SeqIO.SeqRecord):
        logger.debug(f"Mapping DNA record: {record.id}")

        # Get the organism information
        organism_dict = self.map_organism(record)
        organism = Organism.get_or_save(
            taxonomy_id=organism_dict["taxonomy_id"], name=organism_dict["name"]
        )

        dna_infos_dict = self.map_general_infos(record)

        try:
            dna = DNA.get_or_save(
                accession_id=record.id,
                sequence=str(record.seq),
                name=dna_infos_dict["name"],
                seq_length=int(dna_infos_dict["seq_length"]),
                gc_content=float(dna_infos_dict["gc_content"]),
            )
        except Exception as e:
            logger.error(f"Error saving DNA record {record.id}: {e}")
            return

        dna.organism.connect(organism)

        # Here we get the sites informations
        sites = self.map_sites(record)
        self.add_sites(dna, sites)

        # Here we get the regions informations
        regions = self.map_regions(record)
        self.add_regions(dna, regions)

        # check encodings for known proteins
        self.map_cds(dna, record)

    def add_sites(self, dna: DNA, sites: List[dict]):
        for site in sites:
            try:
                site_saving = Site.get_or_save(
                    site_id=site["id"],
                    annotation=Annotation.DNA.value,
                )

                dna.site.connect(
                    site_saving, {"positions": list(range(site["start"], site["end"]))}
                )

            except Exception as e:
                logger.error(
                    f"Error saving site {site['id']} for {dna.accession_id}: {e}"
                )

    def add_regions(self, dna: DNA, regions: List[dict]):
        for region in regions:
            try:
                region_saving = Region.get_or_save(
                    region_id=region["id"],
                    annotation=Annotation.DNA.value,
                )

                dna.region.connect(
                    region_saving, {"start": region["start"], "end": region["end"]}
                )

            except Exception as e:
                logger.error(
                    f"Error saving region {region['id']} for {dna.accession_id}: {e}"
                )

    def map_general_infos(self, seq_record: SeqRecord):
        """
        Extracts general information from a DNA sequence record.
        """
        dna_infos_dict = {}

        dna_infos_dict["name"] = seq_record.name
        dna_infos_dict["seq_length"] = str(len(seq_record.seq))
        dna_infos_dict["gc_content"] = str(self.calculate_gc_content(seq_record.seq))

        return dna_infos_dict

    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculates the GC content of a DNA sequence.
        """
        gc_count = sequence.count("G") + sequence.count("C")
        gc_content = gc_count / len(sequence)
        return gc_content

    def map_organism(self, seq_record: SeqRecord) -> dict:
        """
        Gets the organism name and taxonomy ID from the source data.
        Maps it to an Organism object.
        """

        feature_list = self.get_feature(seq_record, "source")
        if len(feature_list) < 1:
            logger.debug(
                f"Multiple features ({len(feature_list)}) of type `source` found for {seq_record.id}: {feature_list}"
            )
        feature = feature_list[0]

        try:
            if len(feature.qualifiers["db_xref"]) != 1:
                logger.info(
                    f"For {seq_record.id} {feature.qualifiers['db_xref']} taxonomy ID(s) were found, using the first one. Skipping organism assignment"
                )
                return {}

            try:
                taxonomy_id = next(
                    feature
                    for feature in feature.qualifiers["db_xref"]
                    if "taxon" in feature
                )
                if ":" in taxonomy_id:
                    taxonomy_id = taxonomy_id.split(":")[1]
            except StopIteration:
                taxonomy_id = None

        except KeyError:
            logger.debug(f"No taxonomy ID found for {seq_record.id}: {feature}")
            return {}

        try:
            organism_name = feature.qualifiers["organism"]
        except KeyError:
            logger.debug(
                f"No organism name found for {seq_record.id}: {feature_list[0].qualifiers}"
            )
            organism_name = ""

        return {
            "id": taxonomy_id,
            "name": organism_name[0],
            "taxonomy_id": int(taxonomy_id),
        }

    def get_feature(self, seq_record: SeqRecord, feature_type: str) -> list[SeqFeature]:
        """Returns a list of features of a given type from a `Bio.SeqRecord` object."""
        return [
            feature
            for feature in seq_record.features
            if feature.type.lower() == feature_type.lower()
        ]

    def map_sites(self, seq_record: SeqRecord):
        sites_list = []

        sites = self.get_feature(seq_record, "variation")

        for site in sites:
            sites_dict = {}

            try:
                sites_dict["id"] = site.qualifiers["note"][0]
                sites_dict["start"] = int(site.location.start)
                sites_dict["end"] = int(site.location.end)

            except KeyError:
                logger.debug(
                    f"Error mapping site for {seq_record.id}: {site.qualifiers}"
                )

            sites_list.append(sites_dict)

        return sites_list

    def map_regions(self, seq_record: SeqRecord):
        regions_list = []

        regions = self.get_feature(seq_record, "gene")

        for i, region in enumerate(regions):
            try:
                if "db_xref" not in region.qualifiers:
                    db_xref = None
                else:
                    db_xref = region.qualifiers["db_xref"][0]

                regions_list.append(
                    {
                        "id": region.qualifiers["gene"][0],
                        "start": int(region.location.start),
                        "end": int(region.location.end),
                        "cross_reference": db_xref,
                    }
                )

            except KeyError:
                logger.debug(
                    f"Error mapping region for {seq_record.id}: {region.qualifiers}"
                )

        return regions_list

    def map_cds(self, dna, seq_record: SeqRecord):
        cds_list_features = self.get_feature(seq_record, "CDS")

        for cds in cds_list_features:
            try:
                protein_id = cds.qualifiers["protein_id"][0]

                # check if the protein sequence is in database
                protein = Protein.get_or_save(
                    accession_id=protein_id,
                    sequence=cds.qualifiers["translation"][0],
                    seq_length=len(cds.qualifiers["translation"][0]),
                )

                dna.protein.connect(
                    protein,
                    {"start": int(cds.location.start), "end": int(cds.location.end)},
                )

            except KeyError:
                logger.debug(f"Error mapping CDS for {seq_record.id}: {cds.qualifiers}")
