import re
import secrets
from datetime import datetime
from tqdm import tqdm
from typing import List
from Bio import SeqIO, Entrez
from pyEED.core.citation import Citation
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
    accession_id: str, database: str, email: str = None, api_key: str = None
) -> SeqIO.SeqRecord:
    # generate generic mail if none is given
    if email is None:
        email = f"{secrets.token_hex(8)}@gmail.com"

    databases = {"nucleotide", "protein"}
    if database not in databases:
        raise ValueError(f"database must be one of {databases}")

    Entrez.email = email
    Entrez.api_key = api_key

    handle = Entrez.efetch(
        db=database, id=accession_id, rettype="gb", retmode="text", api_key=api_key
    )
    seq_record = SeqIO.read(handle, "genbank")

    handle.close()
    return seq_record


def get_ncbi_entrys(
    accession_ids: List[str], database: str, email: str = None, api_key: str = None
) -> SeqIO.SeqRecord:
    # generate generic mail if none is given
    if email is None:
        email = f"{secrets.token_hex(8)}@gmail.com"

    # Concat accession_ids to string
    if isinstance(accession_ids, list):
        accession_id = ",".join(accession_ids)
    else:
        raise ValueError("'accession_ids' must be a list")

    databases = {"nucleotide", "protein"}
    if database not in databases:
        raise ValueError(f"database must be one of {databases}")

    # Make request
    Entrez.email = email
    Entrez.api_key = api_key

    handle = Entrez.efetch(
        db=database, id=accession_id, rettype="gb", retmode="text", api_key=api_key
    )

    seq_records = []
    for record in tqdm(
        SeqIO.parse(handle, "genbank"),
        desc="⬇️ Fetching protein sequences",
        total=len(accession_ids),
    ):
        seq_records.append(record)

    handle.close()

    return seq_records


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

    try:
        domain = taxonomy[0]
    except IndexError:
        domain = None

    try:
        kingdom = taxonomy[1]
    except IndexError:
        kingdom = None

    try:
        phylum = taxonomy[3]
    except IndexError:
        phylum = None

    try:
        tax_class = taxonomy[5]
    except IndexError:
        tax_class = None

    try:
        order = taxonomy[9]
    except IndexError:
        order = None

    try:
        family = taxonomy[13]
    except IndexError:
        family = None

    try:
        genus = taxonomy[14]
    except IndexError:
        genus = None

    try:
        species = feature.qualifiers["organism"][0].split(" ")[1]
    except IndexError:
        species = None

    return Organism(
        name=feature.qualifiers["organism"][0],
        taxonomy_id=feature.qualifiers["db_xref"][0],
        domain=domain,
        kingdom=kingdom,
        phylum=phylum,
        tax_class=tax_class,
        order=order,
        family=family,
        genus=genus,
        species=species,
    )


def _seqio_to_nucleotide_info(cls, entry: SeqIO):
    """Handel SeqIO entry and return `ProteinSequence`"""

    try:
        sites = []
        regions = []
        for feature in entry.features:
            # TODO: assert that only one protein is in the file
            if feature.type == "Protein":
                if "product" in feature.qualifiers:
                    protein_name = feature.qualifiers["product"][0]
                elif "name" in feature.qualifiers:
                    protein_name = feature.qualifiers["name"][0]
                else:
                    print(f"No name info found in {entry.id}")
                    protein_name = None

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
                    
                if "note" in feature.qualifiers:
                    if feature.qualifiers["note"]:
                        note = feature.qualifiers["note"][0]
                    else:
                        note = ""
                else:
                    note = ""

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
                        note=note,
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
                cds_regions = None

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
    except UnboundLocalError:
        print(f"❌ Sequence {entry.id} was not added, since mapping failed.")


def get_cds_regions(coded_by: dict) -> List[DNARegion]:
    """Extract information about the coding sequence from the 'coded_by' qualifier."""

    cds_pattern = r"\w+\.\d+:\d+\.\.\d+\s?\d+"

    # Extract all regions from the 'coded_by' qualifier
    coded_by = coded_by.replace(">", "")
    coded_by = coded_by.replace("<", "")
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
