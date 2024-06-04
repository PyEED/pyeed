from __future__ import annotations

import json
import logging
from typing import TYPE_CHECKING, List

from pyeed.core import Annotation, Organism

if TYPE_CHECKING:
    from pyeed.core import ProteinRecord

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
LOGGER.addHandler(logging.StreamHandler())


class PDBMapper:
    def __init__(self):
        pass

    def map_pdb_data(self, pdb_data: str) -> List[ProteinRecord]:
        from pyeed.core import ProteinRecord

        records = []

        pdb_entries = json.loads(pdb_data)
        entries = pdb_entries["data"]["entry"]["polymer_entities"]
        try:
            assert len(entries) > 0, f"No entries found in the PDB data {pdb_data}"
        except TypeError:
            return []
        for entry in entries:
            structure_id = entry["rcsb_id"]
            seq_info = entry["rcsb_polymer_entity_container_identifiers"]
            seq_organism = entry["rcsb_entity_source_organism"]
            try:
                seq_id = next(
                    (
                        identifier["database_accession"]
                        for identifier in seq_info["reference_sequence_identifiers"]
                        if "database_accession" in identifier
                    ),
                    None,
                )
            except TypeError:
                continue
            if not seq_id:
                continue

            try:
                tax_id = next(
                    (
                        identifier["ncbi_taxonomy_id"]
                        for identifier in seq_organism
                        if "ncbi_taxonomy_id" in identifier
                    ),
                    None,
                )
            except TypeError:
                tax_id = None

            sequence = entry["entity_poly"]["pdbx_seq_one_letter_code"]

            if tax_id:
                organism = Organism(
                    id=str(tax_id),
                    taxonomy_id=int(tax_id),
                )
            else:
                organism = None

            prot_record = ProteinRecord(
                id=seq_id,
                sequence=sequence,
                organism=organism,
                structure_id=structure_id,
            )

            if entry["rcsb_polymer_entity_feature"]:
                for feature in entry["rcsb_polymer_entity_feature"]:
                    if feature["type"] == "Pfam":
                        for region in feature["feature_positions"]:
                            region = prot_record.add_to_regions(
                                id=feature["feature_id"],
                                name=feature["type"],
                                start=region["beg_seq_id"],
                                end=region["end_seq_id"],
                            )
                            region.add_object_term(Annotation.FAMILY.value)

            try:
                for instance in entry["polymer_entity_instances"]:
                    polymer_instances = instance.get(
                        "rcsb_polymer_instance_feature", None
                    )
                    if not polymer_instances:
                        continue
                    for feature in instance["rcsb_polymer_instance_feature"]:
                        if feature["name"] == "sheet":
                            for strand in feature["feature_positions"]:
                                region = prot_record.add_to_regions(
                                    name=feature["name"],
                                    start=strand["beg_seq_id"],
                                    end=strand["end_seq_id"],
                                )
                                region.add_object_term(Annotation.BETASTRAND.value)
                                region.add_object_term(
                                    "http://edamontology.org/topic_3542"
                                )  # secondary structure

                        elif feature["name"] == "helix":
                            for helix in feature["feature_positions"]:
                                region = prot_record.add_to_regions(
                                    name=feature["name"],
                                    start=helix["beg_seq_id"],
                                    end=helix["end_seq_id"],
                                )
                                region.add_object_term(Annotation.ALPHAHELIX.value)
                                region.add_object_term(
                                    "http://edamontology.org/topic_3542"
                                )  # secondary structure

                        elif feature["name"] == "binding_site":
                            positions = [
                                site["beg_seq_id"]
                                for site in feature["feature_positions"]
                            ]
                            site = prot_record.add_to_sites(
                                name=feature["name"], positions=positions
                            )
                            site.add_object_term(Annotation.BINDING_SITE.value)

            except KeyError:
                pass

            records.append(prot_record)

        return records
