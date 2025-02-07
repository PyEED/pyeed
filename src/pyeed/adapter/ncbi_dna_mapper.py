import io
from typing import Any, List

from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from httpx import Response
from loguru import logger

from pyeed.adapter.primary_db_adapter import PrimaryDBMapper
from pyeed.model import DNA, Annotation, Organism, Protein, Region, Site


class NCBIDNAToPyeed(PrimaryDBMapper):
    """
    Maps DNA sequence entries from NCBI to the PyEED graph object model and saves
    them to the database.
    """

    def add_to_db(self, response: Response) -> None:
        """
        Add a DNA record to the database.

        Args:
            record (SeqIO.SeqRecord): A `Bio.SeqRecord` object.

        Returns:
            None
        """

        with open("tests/data/api_responses/ncbi_dna_QLYQ01000020.txt", "w") as f:
            f.write(response.content.decode())
        assert (
            response.status_code == 200
        ), f"Request to {response.url} failed with status code {response.status_code}"

        records = self.parse_response(response.content)

        for record in records:
            logger.debug(f"Mapping DNA record: {record.id}")

            dna_infos_dict = self.map_general_infos(record)
            if not dna_infos_dict:
                continue

            # Get the organism information
            organism_dict = self.map_organism(record)
            organism = Organism.get_or_save(
                taxonomy_id=organism_dict["taxonomy_id"], name=organism_dict["name"]
            )

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

            dna.organism.connect(organism)

            # Here we get the sites informations
            sites = self.map_sites(record)
            self.add_sites(dna, sites)

            # Here we get the regions informations
            regions = self.map_regions(record)
            self.add_regions(dna, regions)

            # self.map_cds(dna, record)

    def add_sites(self, dna: DNA, sites: List[dict[str, Any]]) -> None:
        """
        Add sites to a DNA record.
        """
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

    def add_regions(self, dna: DNA, regions: List[dict[str, Any]]) -> None:
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

    def map_general_infos(self, seq_record: SeqRecord) -> dict[str, Any]:
        """
        Extracts general information from a DNA sequence record.
        """
        try:
            return dict(
                name=seq_record.name,
                seq_length=len(str(seq_record.seq)),
                gc_content=self.calculate_gc_content(seq_record.seq),
            )
        except UndefinedSequenceError:
            logger.error(f"Undefined sequence error for {seq_record.id}")
            return {}

    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculates the GC content of a DNA sequence.
        """
        gc_count = sequence.count("G") + sequence.count("C")
        gc_content = gc_count / len(sequence)
        return gc_content

    def map_organism(self, seq_record: SeqRecord) -> dict[str, Any]:
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

    def map_sites(self, seq_record: SeqRecord) -> List[dict[str, Any]]:
        sites_list = []

        sites = self.get_feature(seq_record, "variation")

        for site in sites:
            sites_dict = {}

            try:
                sites_dict["id"] = site.qualifiers["note"][0]

                if hasattr(site.location, "start") and hasattr(site.location, "end"):
                    sites_dict["start"] = int(site.location.start)  # type: ignore
                    sites_dict["end"] = int(site.location.end)  # type: ignore
                else:
                    logger.debug(
                        f"Error mapping site for {seq_record.id}: {site.location}"
                    )

            except KeyError:
                logger.debug(
                    f"Error mapping site for {seq_record.id}: {site.qualifiers}"
                )

            sites_list.append(sites_dict)

        return sites_list

    def map_regions(self, seq_record: SeqRecord) -> List[dict[str, Any]]:
        """
        Maps regions to a DNA record.
        """
        regions_list = []

        regions = self.get_feature(seq_record, "gene")

        for region in regions:
            try:
                if "locus_tag" not in region.qualifiers:
                    locus_tag = None
                else:
                    locus_tag = region.qualifiers["locus_tag"][0]

                if hasattr(region.location, "start") and hasattr(
                    region.location, "end"
                ):
                    regions_list.append(
                        {
                            "id": locus_tag,
                            "start": int(region.location.start),  # type: ignore
                            "end": int(region.location.end),  # type: ignore
                            "cross_reference": locus_tag,
                        }
                    )

            except KeyError:
                logger.debug(
                    f"Error mapping region for {seq_record.id}: {region.qualifiers}"
                )

        return regions_list

    def map_cds(self, dna: DNA, seq_record: SeqRecord, protein_id: str) -> None:
        """
        Maps a  CDS to a DNA record.
        """
        # TODO: refactor sequence mapping to protein, allowing to specify the protein id
        cds_list_features = self.get_feature(seq_record, "CDS")

        for cds in cds_list_features:
            if protein_id != cds.qualifiers["protein_id"][0]:
                continue

            # check if the protein sequence is in database
            protein = Protein.get_or_save(
                accession_id=protein_id,
                sequence=cds.qualifiers["translation"][0],
                seq_length=len(cds.qualifiers["translation"][0]),
            )

            if hasattr(cds.location, "start") and hasattr(cds.location, "end"):
                dna.protein.connect(
                    protein,
                    {
                        "start": int(cds.location.start),  # type: ignore
                        "end": int(cds.location.end),  # type: ignore
                    },
                )
            else:
                logger.debug(f"Error mapping CDS for {seq_record.id}: {cds.location}")

    def parse_response(self, content: bytes) -> list[SeqRecord]:
        """Parse NCBI GenBank format response."""
        if not content:
            return []

        try:
            return list(SeqIO.parse(io.StringIO(content.decode()), "gb"))  # type: ignore
        except Exception as e:
            logger.error(f"Failed to parse NCBI DNA response: {e}")
            return []
