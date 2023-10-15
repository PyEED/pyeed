import re
from typing import List, Tuple
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from pyEED.core.nucleotidesequence import NucleotideSequence

from pyEED.core.region import Region
from pyEED.core.site import Site

from ..core.organism import Organism


def SeqIO_to_pyeed(entry: SeqIO):
    """Handel SeqIO entry and return pyeed object."""

    if entry.annotations["molecule_type"] == "protein":
        return _seqio_to_protein_sequence(entry)

    elif entry.annotations["molecule_type"] == "DNA":
        raise NotImplementedError("DNA is not implemented yet.")

    else:
        raise ValueError(
            f"{entry.id} of type {entry.annotations['molecule_type']} is not 'protein' or 'DNA'."
        )


def _seqio_to_protein_sequence(cls, entry: SeqIO):
    """Handel SeqIO entry and return `ProteinSequence`"""

    if "db_source" in entry.annotations:
        if "pdb" in entry.annotations["db_source"]:
            pdb_id = entry.id
        else:
            pdb_id = None

    sites = []
    protein_regions = []
    for feature in entry.features:
        # TODO: assert that only one protein is in the file
        if feature.type == "Protein":
            if "product" in feature.qualifiers:
                protein_name = feature.qualifiers["product"][0]

            if "calculated_mol_wt" in feature.qualifiers:
                mol_weight = feature.qualifiers["calculated_mol_wt"][0]
            else:
                mol_weight = None

            if "EC_number" in feature.qualifiers:
                ec_number = feature.qualifiers["EC_number"][0]
            else:
                ec_number = None

        if feature.type == "source":
            organism = Organism(
                name=feature.qualifiers["organism"][0],
                taxonomy_id=feature.qualifiers["db_xref"][0].split(":")[1],
                id=feature.qualifiers["db_xref"][0],
            )

        if feature.type == "Region":
            protein_regions.append(
                Region(
                    name=feature.qualifiers["region_name"][0],
                    start=int(feature.location.start),
                    end=int(feature.location.end),
                    cross_reference=feature.qualifiers["db_xref"][0],
                    note=feature.qualifiers["note"][0],
                )
            )

        if feature.type == "Site":
            site_type = feature.qualifiers["site_type"][0]

            if "note" in feature.qualifiers:
                name = feature.qualifiers["note"][0]
            else:
                name = site_type

            sites.append(
                Site(
                    name=name,
                    positions=[loc for loc in feature.location],
                    cross_reference=feature.qualifiers["db_xref"][0],
                    type=site_type,
                )
            )

        if feature.type == "CDS":
            cds_id, cds_regions = get_cds_info(feature.qualifiers["coded_by"][0])

            # TODO: find better way to discriminate between 'gene' and 'locus_tag'
            if "gene" in feature.qualifiers:
                gene_id = feature.qualifiers["gene"][0]
            else:
                try:
                    gene_id = feature.qualifiers["locus_tag"][0]
                except KeyError:
                    gene_id = None

            coding_sequence = NucleotideSequence(
                id=cds_id, regions=cds_regions, gene_id=gene_id
            )

    return cls(
        id=entry.id,
        name=protein_name,
        sequence=str(entry.seq),
        ec_number=ec_number,
        mol_weight=mol_weight,
        organism=organism,
        sites=sites,
        regions=protein_regions,
        pdb_id=pdb_id,
        coding_sequence=coding_sequence,
    )


def get_cds_info(coded_by: dict) -> Tuple[str, List[Region]]:
    """Extract information about the coding sequence from the 'coded_by' qualifier."""

    cds_pattern = r"\w+\.\d+:\d+\.\.\d+\s?\d+"

    # Extract all regions from the 'coded_by' qualifier
    cds_regions = re.findall(cds_pattern, coded_by)
    cds_regions = [region.replace(" ", "") for region in cds_regions]

    # Extract the reference id from the first region
    reference_ids = [region.split(":")[0] for region in cds_regions]
    if not all([reference_id == reference_ids[0] for reference_id in reference_ids]):
        print("nucleotide sequence references are not identical.")

    # Extract the start and end position of each region
    regions = []
    cds_ranges = [region.split(":")[1] for region in cds_regions]
    for region in cds_ranges:
        start, end = region.split("..")
        regions.append(Region(id=reference_ids[0], start=int(start), end=int(end)))

    return reference_ids[0], regions


def extract_nucleotide_seq(entry: SeqIO, nucleotide_sequence: NucleotideSequence):
    """Handel nucleotide SeqIO entry and map it to `NucleotideSequence`"""

    feature_regions = set()
    for region in nucleotide_sequence.regions:
        feature_regions.add(region.start)
        feature_regions.add(region.end)

    for feature in entry.features:
        if feature.type == "source":
            nucleotide_sequence.molecule_type = feature.qualifiers["mol_type"][0]

        if feature.type == "CDS":
            if isinstance(feature.location, CompoundLocation):
                locations = set()
                parts = feature.location.parts
                for part in parts:
                    # TODO: investigate reason for +1
                    locations.add(int(part.start) + 1)
                    locations.add(int(part.end))

            if isinstance(feature.location, FeatureLocation):
                locations = {int(feature.location.start) + 1, int(feature.location.end)}

            if feature_regions == locations:
                nucleotide_sequence.sequence = str(feature.location.extract(entry.seq))
