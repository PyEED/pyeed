from __future__ import annotations

import io
import logging
import re
from typing import TYPE_CHECKING, List

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from pydantic import ValidationError
from pyeed.core.annotation import Annotation
from pyeed.core.organism import Organism
from pyeed.core.region import Region

if TYPE_CHECKING:
    from pyeed.core.dnarecord import DNARecord


LOGGER = logging.getLogger(__name__)


class NCBIDNAMapper:
    def __init__(self):
        pass

    def _to_seq_records(self, responses: List[str]) -> List[SeqRecord]:
        """
        Converts the fetched data to a list of `Bio.SeqRecord.SeqRecord` objects.
        """
        records = []
        for response in responses:
            LOGGER.debug(f"Converting response to SeqRecord: {response[0]}")
            records.extend(SeqIO.parse(io.StringIO(response[0]), "gb"))

        return records

    def map(self, responses: List[str]) -> List[DNARecord]:
        """
        Maps the fetched data to an instance of the `DNARecord` class.
        """

        from pyeed.core.dnarecord import DNARecord

        seq_records = self._to_seq_records(responses)

        dna_infos = []
        for record in seq_records:
            dna_info = DNARecord(id=record.id, sequence=str(record.seq))

            try:
                dna_info.organism = Organism(**self.map_organism(record))
            except ValidationError as e:
                LOGGER.error(
                    f"Error mapping organism for {record.id}: {e.errors()} {e.json()}"
                )
                continue

            dna_info = self.map_general_infos(record, dna_info)

            dna_info = self.map_regions(record, dna_info)

            dna_info = self.map_sites(record, dna_info)

            # dna_info = self.map_cds(record, dna_info)

            dna_infos.append(dna_info)

        return dna_infos

    def map_organism(self, seq_record: SeqRecord) -> dict:
        """
        Gets the organism name and taxonomy ID from the source data.
        Maps it to an Organism object.
        """

        feature = self.get_feature(seq_record, "source")
        if len(feature) < 1:
            LOGGER.debug(
                f"Multiple features ({len(feature)}) of type `source` found for {seq_record.id}: {feature}"
            )
        feature = feature[0]

        try:
            if len(feature.qualifiers["db_xref"]) != 1:
                LOGGER.info(
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
            LOGGER.debug(f"No taxonomy ID found for {seq_record.id}: {feature}")
            return {}

        try:
            organism_name = feature.qualifiers["organism"]
        except KeyError:
            LOGGER.debug(
                f"No organism name found for {seq_record.id}: {feature[0].qualifiers}"
            )
            organism_name = ""

        return {"id": taxonomy_id, "name": organism_name[0], "taxonomy_id": taxonomy_id}

    def map_general_infos(self, seq_record: SeqRecord, dna_info: DNARecord):
        """
        Maps general information from a `Bio.SeqRecord` to a `DNARecord` object.
        """

        dna_info.name = seq_record.name
        dna_info.seq_length = len(seq_record.seq)
        # GC-content=(A+T+G+C)(G+C)​×100
        dna_info.gc_content = (
            (seq_record.seq.count("G") + seq_record.seq.count("C"))
            / dna_info.seq_length
        ) * 100

        return dna_info

    def map_sites(self, seq_record: SeqRecord, dna_info: DNARecord):
        """
        Maps site data from a `Bio.SeqRecord` to a `DNARecord` object.
        """

        sites = self.get_feature(seq_record, "variation")
        for site in sites:
            try:
                dna_info.sites.append(
                    Region(
                        id=site.qualifiers["note"][0],
                        start=int(site.location.start),
                        end=int(site.location.end),
                    )
                )
            except KeyError:
                LOGGER.debug(
                    f"Error mapping site for {seq_record.id}: {site.qualifiers}"
                )

        return dna_info

    def map_regions(self, seq_record: SeqRecord, dna_info: DNARecord):
        """Maps region data from a `Bio.SeqRecord` to a `DNARecord` object."""

        regions = self.get_feature(seq_record, "gene")
        for region in regions:
            try:
                if "db_xref" not in region.qualifiers:
                    db_xref = None
                else:
                    db_xref = region.qualifiers["db_xref"][0]

                dna_info.regions.append(
                    Region(
                        id=region.qualifiers["gene"][0],
                        start=int(region.location.start),
                        end=int(region.location.end),
                        cross_reference=db_xref,
                    )
                )

            except KeyError:
                LOGGER.debug(
                    f"Error mapping region for {seq_record.id}: {region.qualifiers}"
                )

        return dna_info

    def map_cds(self, seq_record: SeqRecord, protein_record: ProteinRecord):
        """Maps coding sequence data from a `Bio.SeqRecord` to a `ProteinRecord` object."""

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

    def get_feature(self, seq_record: SeqRecord, feature_type: str) -> List[SeqFeature]:
        """Returns a list of features of a given type from a `Bio.SeqRecord` object."""
        return [
            feature
            for feature in seq_record.features
            if feature.type.lower() == feature_type.lower()
        ]
