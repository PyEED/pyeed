import io
import re
from typing import Any, Dict, Generator, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from httpx import Response
from loguru import logger

from pyeed.adapter.primary_db_adapter import PrimaryDBMapper
from pyeed.model import Annotation, Organism, Protein, Region, Site


class NCBIProteinToPyeed(PrimaryDBMapper):
    def get_feature(self, seq_record: SeqRecord, feature_type: str) -> List[SeqFeature]:
        """Extract features of a given type from a SeqRecord."""
        features = [
            feature
            for feature in seq_record.features
            if feature.type.lower() == feature_type.lower()
        ]
        if feature_type.lower() == "source" and len(features) != 1:
            logger.warning(
                f"Record {seq_record.id}: Expected 1 'source' feature, found {len(features)}."
            )
        return features

    def map_organism(self, seq_record: SeqRecord) -> Tuple[Any, Any]:
        """
        Extract organism name and taxonomy ID from the source feature.
        """
        source_features = self.get_feature(seq_record, "source")
        if not source_features:
            logger.warning(f"Record {seq_record.id}: No 'source' feature found.")
            return (None, None)
        feature = source_features[0]

        taxonomy_id = None
        try:
            db_xrefs = feature.qualifiers.get("db_xref", [])
            if len(db_xrefs) != 1:
                logger.debug(
                    f"Record {seq_record.id}: Multiple taxonomy references found; using first valid one."
                )
            for ref in db_xrefs:
                if ref.startswith("taxon:"):
                    taxonomy_id = int(ref.split(":")[1])
                    logger.debug(f"Record {seq_record.id}: Taxonomy ID: {taxonomy_id}")
                    break
            if taxonomy_id is None:
                logger.debug(
                    f"Record {seq_record.id}: No valid taxonomy ID found in db_xref {db_xrefs}."
                )
                return (None, None)
        except KeyError as e:
            logger.error(
                f"Record {seq_record.id}: Exception extracting taxonomy ID: {e}."
            )
            return (None, None)

        try:
            organism_names = feature.qualifiers.get("organism", [])
            if not organism_names:
                logger.warning(
                    f"Record {seq_record.id}: 'organism' qualifier missing in source feature."
                )
                organism_name = None
            else:
                organism_name = organism_names[0]
        except Exception as e:
            logger.error(
                f"Record {seq_record.id}: Exception extracting organism name: {e}."
            )
            organism_name = None

        logger.debug(
            f"Record {seq_record.id}: Mapped organism '{organism_name}' with taxonomy ID {taxonomy_id}."
        )
        return organism_name, taxonomy_id

    def map_protein(self, seq_record: SeqRecord) -> Dict[str, Any]:
        """Map protein information from a SeqRecord to a dictionary."""
        protein_info: Dict[str, Any] = {}
        protein_features = self.get_feature(seq_record, "Protein")
        if not protein_features:
            logger.warning(f"Record {seq_record.id}: No protein feature found.")
            return protein_info

        if len(protein_features) > 1:
            logger.warning(
                f"Record {seq_record.id}: Multiple protein features found; using the first."
            )

        protein = protein_features[0]
        try:
            protein_info["name"] = protein.qualifiers.get("product", [None])[0]
            if protein_info["name"] is None:
                raise KeyError
        except KeyError:
            logger.debug(
                f"Record {seq_record.id}: 'product' qualifier missing; trying 'name' qualifier."
            )
            protein_info["name"] = protein.qualifiers.get("name", [None])[0]
            if protein_info["name"] is None:
                logger.debug(
                    f"Record {seq_record.id}: Protein name not provided; setting to None."
                )

        try:
            protein_info["mol_weight"] = float(
                protein.qualifiers.get("calculated_mol_wt", [None])[0]
            )
        except (KeyError, ValueError, TypeError):
            logger.debug(
                f"Record {seq_record.id}: Molecular weight missing or invalid; setting to None."
            )
            protein_info["mol_weight"] = None

        try:
            protein_info["ec_number"] = protein.qualifiers.get("EC_number", [None])[0]
        except KeyError:
            logger.warning(
                f"Record {seq_record.id}: EC number not provided; setting to None."
            )
            protein_info["ec_number"] = None

        logger.debug(
            f"Record {seq_record.id}: Mapped protein with name '{protein_info.get('name')}' to 'Protein' object."
        )
        return protein_info

    def map_sites(self, seq_record: SeqRecord) -> List[Dict[str, Any]]:
        """Map site features from a SeqRecord into a list of dictionaries."""
        sites_list = []
        site_features = self.get_feature(seq_record, "site")
        for site in site_features:
            try:
                sites_dict = {
                    "type": site.qualifiers["site_type"][0],
                    "positions": [int(part.start) for part in site.location.parts],  # type: ignore
                    "cross_ref": site.qualifiers["db_xref"][0],
                }
                sites_list.append(sites_dict)
            except KeyError:
                logger.warning(
                    f"Record {seq_record.id}: Incomplete site data; skipping site."
                )
        logger.debug(f"Record {seq_record.id}: Mapped {len(sites_list)} site(s).")
        return sites_list

    def add_sites(self, sites_list: List[Dict[str, Any]], protein: Protein) -> None:
        """Connect mapped sites to a protein."""
        for site_dict in sites_list:
            site = Site.get_or_save(
                name=site_dict["type"],
                annotation=Annotation.PROTEIN.value,
            )
            protein.site.connect(site, {"positions": site_dict["positions"]})
        logger.debug(
            f"Connected {len(sites_list)} site(s) to protein {protein.accession_id}."
        )

    def map_cds(self, seq_record: SeqRecord) -> Optional[List[Dict[str, Any]]]:
        """Map CDS features from a SeqRecord."""
        cds_features = self.get_feature(seq_record, "CDS")
        if len(cds_features) > 1:
            logger.debug(
                f"Record {seq_record.id}: Multiple CDS features found; using the first."
            )
        try:
            cds_feature = cds_features[0]
        except IndexError:
            logger.debug(f"Record {seq_record.id}: No CDS feature found.")
            return None

        if "coded_by" not in cds_feature.qualifiers:
            logger.debug(
                f"Record {seq_record.id}: CDS feature missing 'coded_by' qualifier; skipping CDS mapping."
            )
            return None

        logger.debug(f"Record {seq_record.id}: Processing CDS feature with qualifiers.")
        return self.get_cds_regions(cds_feature)

    @staticmethod
    def get_cds_regions(seq_feature: SeqFeature) -> List[Dict[str, Any]]:
        """Extract CDS regions from a CDS feature."""
        coded_by = seq_feature.qualifiers["coded_by"][0]
        cds_pattern = r"\w+\.\d+:\d+\.\.\d+\s?\d+"
        coded_by = coded_by.replace(">", "").replace("<", "")
        cds_regions = re.findall(cds_pattern, coded_by)
        cds_regions = [region.replace(" ", "") for region in cds_regions]
        reference_ids = [region.split(":")[0] for region in cds_regions]
        if not all(ref == reference_ids[0] for ref in reference_ids):
            logger.warning(
                f"Inconsistent nucleotide references found: {reference_ids}."
            )
        cds_ranges = [region.split(":")[1] for region in cds_regions]
        regions = []
        for region_str in cds_ranges:
            start, end = region_str.split("..")
            regions.append(
                {
                    "start": int(start),
                    "end": int(end),
                    "type": Annotation.CODING_SEQ.value,
                    "id": reference_ids[0],
                }
            )
        logger.debug(f"Extracted {len(regions)} CDS region(s).")
        return regions

    def map_regions(self, seq_record: SeqRecord) -> List[Dict[str, Any]]:
        """Map region features from a SeqRecord."""
        region_features = self.get_feature(seq_record, "region")
        regions_list = []
        for region in region_features:
            try:
                db_xref = region.qualifiers.get("db_xref", [None])[0]
                regions_list.append(
                    {
                        "id": region.qualifiers["region_name"][0],
                        "start": int(region.location.start),  # type: ignore
                        "end": int(region.location.end),  # type: ignore
                        "type": region.qualifiers["note"][0],
                        "cross_reference": db_xref,
                    }
                )
            except KeyError:
                logger.warning(
                    f"Record {seq_record.id}: Incomplete region data; skipping region."
                )
        logger.debug(f"Record {seq_record.id}: Mapped {len(regions_list)} region(s).")
        return regions_list

    def add_regions(self, regions_list: List[Dict[str, Any]], protein: Protein) -> None:
        """Connect mapped regions to a protein."""
        for region_dict in regions_list:
            region = Region.get_or_save(
                region_id=region_dict["id"],
                annotation=Annotation.PROTEIN.value,
            )
            protein.region.connect(
                region, {"start": region_dict["start"], "end": region_dict["end"]}
            )
        logger.debug(
            f"Connected {len(regions_list)} region(s) to protein {protein.accession_id}."
        )

    def add_to_db(self, response: Response) -> None:
        """Process the response, map record data, and add/update the protein in the database."""
        records = list(self.parse_response(response))
        if not records:
            logger.error("Response parsing failed: no records found.")
            return

        for record in records:
            logger.debug(f"Processing NCBI protein record {record.id}")

            organism_name, taxonomy_id = self.map_organism(record)
            organism = Organism.get_or_save(taxonomy_id=taxonomy_id, name=organism_name)

            protein_info_dict = self.map_protein(record)
            protein_data = {
                "sequence": str(record.seq),
                "mol_weight": protein_info_dict.get("mol_weight"),
                "ec_number": protein_info_dict.get("ec_number"),
                "name": protein_info_dict.get("name"),
                "seq_length": len(record.seq),
                "accession_id": str(record.id),
            }

            # Create or update protein
            protein = Protein.nodes.get_or_none(accession_id=record.id)
            if protein:
                for key, value in protein_data.items():
                    setattr(protein, key, value)
                protein.save()
            else:
                protein = Protein(**protein_data)
                protein.save()

            protein.organism.connect(organism)

            # Add features
            sites_list = self.map_sites(record)
            self.add_sites(sites_list, protein)

            cds_regions = self.map_cds(record)
            if cds_regions:
                for region in cds_regions:
                    protein.nucleotide_id = region["id"]
                    protein.nucleotide_start = region["start"]
                    protein.nucleotide_end = region["end"]
                    protein.save()

            regions_list = self.map_regions(record)
            self.add_regions(regions_list, protein)

            logger.info(f"Added/updated NCBI protein {record.id} in database")

    def parse_response(self, response: Response) -> Generator[SeqRecord, None, None]:
        """Parse a GenBank-format response from NCBI."""
        if not response.content:
            logger.error("Empty response received from NCBI.")
            return []  # type: ignore
        try:
            return SeqIO.parse(io.StringIO(response.content.decode()), "gb")  # type: ignore
        except Exception as e:
            logger.error(f"Failed to parse GenBank response: {e}")
            return []  # type: ignore
