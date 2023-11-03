import re
import secrets
from typing import List
from Bio import SeqIO, Entrez
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from pyEED.core.dnaregion import DNARegion
from pyEED.core.dnaregiontype import DNARegionType
from pyEED.core.proteinregion import ProteinRegion
from pyEED.core.proteinregiontype import ProteinRegionType
from pyEED.core.proteinsitetype import ProteinSiteType

from pyEED.core.site import Site
from pyEED.core.span import Span

from ..core.organism import Organism

GENEID_PATTERN = r"GeneID:\d+"


def get_ncbi_entry(
    accession_id: str, database: str, email: str = None
) -> SeqIO.SeqRecord:
    # generate generic mail if none is given
    if email is None:
        email = f"{secrets.token_hex(8)}@gmail.com"

    Entrez.email = email

    databases = {"nucleotide", "protein"}
    if database not in databases:
        raise ValueError(f"database must be one of {databases}")

    handle = Entrez.efetch(db=database, id=accession_id, rettype="gb", retmode="text")
    seq_record = SeqIO.read(handle, "genbank")

    handle.close()
    return seq_record


def SeqIO_to_pyeed(entry: SeqIO):
    """Handel SeqIO entry and return pyeed object."""

    if entry.annotations["molecule_type"] == "protein":
        return _seqio_to_nucleotide_info(entry)

    elif entry.annotations["molecule_type"] == "DNA":
        raise NotImplementedError("DNA is not implemented yet.")

    else:
        raise ValueError(
            f"{entry.id} of type {entry.annotations['molecule_type']} is not 'protein' or 'DNA'."
        )


def _seqio_to_dna_info(cls, entry: SeqIO):
    regions = []
    for feature in entry.features:
        if feature.type == "gene":
            regions.append(
                DNARegion(
                    name=feature.qualifiers["gene"][0],
                    spans=[
                        Span(
                            start=int(feature.location.start),
                            end=int(feature.location.end),
                        )
                    ],
                    cross_reference=re.findall(
                        GENEID_PATTERN, "".join(feature.qualifiers["db_xref"])
                    )[0],
                    note=feature.qualifiers["note"][0],
                    type=DNARegionType.GENE,
                )
            )

        if feature.type == "CDS":
            regions.append(
                DNARegion(
                    name=feature.qualifiers["product"][0],
                    spans=[
                        Span(
                            start=int(feature.location.start),
                            end=int(feature.location.end),
                        )
                    ],
                    cross_reference=feature.qualifiers["protein_id"][0],
                    note=feature.qualifiers["note"][0],
                    type=DNARegionType.CODING_SEQUENCE,
                )
            )

        if feature.type == "source":
            organism = get_organism(entry.annotations, feature)

    return cls(
        source_id=entry.id,
        name=entry.description,
        sequence=str(entry.seq),
        regions=regions,
        organism=organism,
    )


def get_organism(annotations, feature) -> Organism:
    taxonomy = annotations["taxonomy"]

    return Organism(
        name=feature.qualifiers["organism"][0],
        taxonomy_id=feature.qualifiers["db_xref"][0],
        domain=taxonomy[0],
        kingdom=taxonomy[1],
        phylum=taxonomy[3],
        tax_class=taxonomy[5],
        order=taxonomy[9],
        family=taxonomy[13],
        genus=taxonomy[14],
    )


def _seqio_to_nucleotide_info(cls, entry: SeqIO):
    """Handel SeqIO entry and return `ProteinSequence`"""

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
            organism = get_organism(entry.annotations, feature)

        if feature.type == "Region":
            if "db_xref" in feature.qualifiers:
                cross_reference = feature.qualifiers["db_xref"][0]
            else:
                cross_reference = None

            regions.append(
                ProteinRegion(
                    name=feature.qualifiers["region_name"][0],
                    spans=[
                        Span(
                            start=int(feature.location.start),
                            end=int(feature.location.end),
                        )
                    ],
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
                    cross_ref=cross_reference,
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
        coding_sequence_ref=cds_regions,
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
    cds_ranges = [region.split(":")[1] for region in cds_regions]
    for region in cds_ranges:
        start, end = region.split("..")
        span = Span(start=int(start), end=int(end))

        region = DNARegion(
            id=reference_ids[0],
            spans=[span],
            type=DNARegionType.CODING_SEQUENCE,
        )

    return region
