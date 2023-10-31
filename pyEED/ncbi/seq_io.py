import re
from typing import List, Tuple
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from pyEED.core.dnainfo import DNAInfo
from pyEED.core.dnaregion import DNARegion
from pyEED.core.dnaregiontype import DNARegionType
from pyEED.core.proteinregion import ProteinRegion
from pyEED.core.proteinregiontype import ProteinRegionType
from pyEED.core.proteinsitetype import ProteinSiteType

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
    else:
        pdb_id = None

    sites = []
    regions = []
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
                taxonomy_id=feature.qualifiers["db_xref"][0],
            )

        if feature.type == "Region":
            if "db_xref" in feature.qualifiers:
                cross_reference = feature.qualifiers["db_xref"][0]
            else:
                cross_reference = None

            regions.append(
                ProteinRegion(
                    name=feature.qualifiers["region_name"][0],
                    start=int(feature.location.start),
                    end=int(feature.location.end),
                    cross_reference=cross_reference,
                    note=feature.qualifiers["note"][0],
                )
            )

        if feature.type == "Site":
            site_type = feature.qualifiers["site_type"][0].lower()

            if "note" in feature.qualifiers:
                name = feature.qualifiers["note"][0]
            else:
                name = site_type

            if "db_xref" in feature.qualifiers:
                cross_reference = feature.qualifiers["db_xref"][0]
            else:
                cross_reference = None

            sites.append(
                Site(
                    name=name,
                    positions=[loc for loc in feature.location],
                    cross_reference=cross_reference,
                    type=ProteinSiteType.match_string(site_type),
                )
            )

        if feature.type == "CDS":
            cds_regions = get_cds_regions(feature.qualifiers["coded_by"][0])

        if "CDS" not in [feature.type for feature in entry.features]:
            cds_regions = []

        if "Protein" not in [feature.type for feature in entry.features]:
            protein_name = entry.description
            ec_number = None
            mol_weight = None

    return cls(
        source_id=entry.id,
        name=protein_name,
        sequence=str(entry.seq),
        ec_number=ec_number,
        mol_weight=mol_weight,
        organism=organism,
        sites=sites,
        regions=regions,
        cds_references=cds_regions,
    )


def get_cds_regions(coded_by: dict) -> List[DNARegion]:
    """Extract information about the coding sequence from the 'coded_by' qualifier."""

    cds_pattern = r"\w+\.\d+:\d+\.\.\d+\s?\d+"

    # Extract all regions from the 'coded_by' qualifier
    coded_by = coded_by.replace(">", "")
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
        regions.append(
            DNARegion(
                id=reference_ids[0],
                start=int(start),
                end=int(end),
                type=DNARegionType.CODING_SEQUENCE,
            )
        )

    return regions


# def extract_nucleotide_seq(entry: SeqIO, nucleotide_sequence: NucleotideSequence):
#     """Handel nucleotide SeqIO entry and map it to `NucleotideSequence`"""

#     feature_regions = set()
#     for region in nucleotide_sequence.regions:
#         feature_regions.add(region.start)
#         feature_regions.add(region.end)

#     for feature in entry.features:
#         if feature.type == "source":
#             nucleotide_sequence.molecule_type = feature.qualifiers["mol_type"][0]

#         if feature.type == "CDS":
#             if isinstance(feature.location, CompoundLocation):
#                 locations = set()
#                 parts = feature.location.parts
#                 for part in parts:
#                     # TODO: investigate reason for +1
#                     locations.add(int(part.start) + 1)
#                     locations.add(int(part.end))

#             if isinstance(feature.location, FeatureLocation):
#                 locations = {int(feature.location.start) + 1, int(feature.location.end)}

#             if feature_regions == locations:
#                 nucleotide_sequence.sequence = str(feature.location.extract(entry.seq))
