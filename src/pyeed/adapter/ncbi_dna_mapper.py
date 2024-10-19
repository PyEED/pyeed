import io
import re
from abc import abstractmethod
from collections import defaultdict
from typing import Generic, TypeVar, List
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from loguru import logger

from pyeed.model import Annotation, GOAnnotation, Organism, Protein, Site, Region, DNA

T = TypeVar("T")

class PrimaryDBtoPyeed(Generic[T]):
    @abstractmethod
    def add_to_db(self, data: dict):
        pass

class NCBIDNAToPyeed(PrimaryDBtoPyeed[DNA]):

    def __init__(self):
        pass

    def add_to_db(self, record: SeqIO.SeqRecord):

        logger.info(f"Mapping DNA record: {record.id}")

        # Here we get the organism information
        organism_dict = self.map_organism(record)
        organism = Organism.get_or_save(
            taxonomy_id=organism_dict["taxonomy_id"], name=organism_dict["name"]
        )

        # Here we get the nucleotide sequence informations (general informations)
        dna_infos_dict = self.map_general_infos(record)


        try:
            dna = DNA.get_or_save(
                accession_id=record.id,
                sequence=str(record.seq),
                name=dna_infos_dict["name"],
                seq_length=dna_infos_dict["seq_length"],
                gc_content=dna_infos_dict["gc_content"],
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

    def add_sites(self, dna: DNA, sites: List[dict]):
        for site in sites:
            try:
                site = Site.get_or_save(
                    site_id=site["id"],
                    # DANGERDANGER
                    annotation=Annotation.ACTIVE_SITE.value,
                )

                dna.sites.connect(site, {'position': list(range(site["start"], site["end"]))})

            except Exception as e:
                logger.error(f"Error saving site {site['id']} for {dna.accession_id}: {e}")

    def add_regions(self, dna: DNA, regions: List[dict]):

        for region in regions:
            try:
                region_saving = Region.get_or_save(
                    region_id=region["id"],
                    # DANGERDANGER
                    annotation=Annotation.ACTIVE_SITE.value,
                )

                dna.region.connect(region_saving, {"start": region["start"], "end": region["end"]})

            except Exception as e:
                logger.error(f"Error saving region {region['id']} for {dna.accession_id}: {e}")

    def map_general_infos(self, seq_record: SeqRecord):

        dna_infos_dict = {}

        dna_infos_dict['name'] = seq_record.name
        dna_infos_dict['seq_length'] = len(seq_record.seq)
        # GC-content=(A+T+G+C)(G+C)​×100
        dna_infos_dict['gc_content'] = (
            (seq_record.seq.count("G") + seq_record.seq.count("C"))
            / dna_infos_dict['seq_length']
        ) * 100

        return dna_infos_dict

    def map_organism(self, seq_record: SeqRecord) -> dict:
        """
        Gets the organism name and taxonomy ID from the source data.
        Maps it to an Organism object.
        """

        feature = self.get_feature(seq_record, "source")
        if len(feature) < 1:
            logger.debug(
                f"Multiple features ({len(feature)}) of type `source` found for {seq_record.id}: {feature}"
            )
        feature = feature[0]

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
                f"No organism name found for {seq_record.id}: {feature[0].qualifiers}"
            )
            organism_name = ""

        return {"id": taxonomy_id, "name": organism_name[0], "taxonomy_id": int(taxonomy_id)}
        
    def get_feature(self, seq_record: SeqRecord, feature_type: str) -> SeqFeature:
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

        for region in regions:
            try:
                if "db_xref" not in region.qualifiers:
                    db_xref = None
                else:
                    db_xref = region.qualifiers["db_xref"][0]

                regions_list.append({
                    "id": region.qualifiers["gene"][0],
                    "start": int(region.location.start),
                    "end": int(region.location.end),
                    "cross_reference": db_xref
                })

            except KeyError:
                logger.debug(
                    f"Error mapping region for {seq_record.id}: {region.qualifiers}"
                )

        return regions_list