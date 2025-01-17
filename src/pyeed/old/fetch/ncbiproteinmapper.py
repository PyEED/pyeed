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
from pyeed.core.dnarecord import DNARecord
from pyeed.core.organism import Organism
from pyeed.core.region import Region

if TYPE_CHECKING:
    from pyeed.core.proteinrecord import ProteinRecord

LOGGER = logging.getLogger(__name__)


class NCBIProteinMapper:
    def __init__(self):
        pass

    def _to_seq_records(self, responses: List[str]) -> List[SeqRecord]:
        """
        Converts the fetched data to a list of `Bio.SeqRecord.SeqRecord` objects.
        """
        records = []
        for response in responses:
            records.extend(SeqIO.parse(io.StringIO(response), "gb"))

        return records

    def map(self, responses: List[str]) -> List[ProteinRecord]:
        """
        Maps the fetched data to an instance of the `ProteinInfo` class.
        """

        from pyeed.core.proteinrecord import ProteinRecord

        seq_records = self._to_seq_records(responses)

        protein_infos = []
        for record in seq_records:
            protein_info = ProteinRecord(id=record.id, sequence=str(record.seq))

            try:
                protein_info.organism = Organism(**self.map_organism(record))
            except ValidationError as e:
                LOGGER.error(
                    f"Error mapping organism for {record.id}: {e.errors()} {e.json()}"
                )
                continue

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

    def map_protein(self, seq_record: SeqRecord, protein_info: ProteinRecord):
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

    def map_regions(self, seq_record: SeqRecord, protein_info: ProteinRecord):
        """Maps region data from a `Bio.SeqRecord` to a `ProteinInfo` object."""

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
                LOGGER.debug(
                    f"Incomplete region data found for {seq_record.id}: {region.qualifiers}, skipping region"
                )

        return protein_info

    def map_sites(self, seq_record: SeqRecord, protein_info: ProteinRecord):
        """Maps site data from a `Bio.SeqRecord` to a `ProteinInfo` object."""

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
