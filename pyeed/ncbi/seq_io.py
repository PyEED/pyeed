import logging
import logging.config
import secrets
from pathlib import Path
from tqdm import tqdm
from typing import List
from Bio import SeqIO, Entrez


GENEID_PATTERN = r"GeneID:\d+"

path_config = Path(__file__).parent.parent.parent / "logging.conf"
logger = logging.getLogger("pyeed")


def get_ncbi_entry(
    accession_id: str,
    database: str,
    retmode: str,
    email: str = None,
    api_key: str = None,
) -> SeqIO.SeqRecord:

    # generate generic mail if none is given
    if email is None:
        email = f"{secrets.token_hex(8)}@gmail.com"

    databases = {"nucleotide", "protein", "taxonomy"}
    if database not in databases:
        raise ValueError(f"Database must be one of {databases}")

    Entrez.email = email
    Entrez.api_key = api_key

    logger.debug(
        f"Fetching {accession_id} from {database} with retmode {retmode} email {email} api_key {api_key}"
    )

    with Entrez.efetch(
        db=database, id=accession_id, rettype="gb", retmode=retmode, api_key=api_key
    ) as handle:
        seq_record = SeqIO.read(handle, "genbank")

    return seq_record


def get_ncbi_entrys(
    accession_ids: List[str],
    database: str,
    retmode: str,
    email: str = None,
    api_key: str = None,
) -> SeqIO.SeqRecord:
    # generate generic mail if none is given
    if email is None:
        email = f"{secrets.token_hex(8)}@gmail.com"

    # database handle
    databases = {"nucleotide", "protein"}
    if database not in databases:
        raise ValueError(f"Database must be one of {databases}")

    # check if accession_ids is a list
    if not isinstance(accession_ids, list):
        raise ValueError("'accession_ids' must be a list")

    # dissect the list in chunks of 100
    accession_sets = [
        accession_ids[i : i + 100] for i in range(0, len(accession_ids), 100)
    ]

    logger.debug(
        f"Subsetted a total of {len(accession_ids)} accession_ids into {len(accession_sets)} sets."
        + f" The last set contains {len(accession_sets[-1])} entries"
    )

    seq_records = []
    with tqdm(
        total=len(accession_ids),
        desc="⬇️ Fetching protein sequences",
    ) as pbar:
        for subset_count, accession_set in enumerate(accession_sets):
            accession_id = ",".join(accession_set)

            # Make request
            Entrez.email = email
            Entrez.api_key = api_key

            logger.debug(
                f"Fetching subset {subset_count+1} of {len(accession_sets)} with {len(accession_set)} entries."
            )

            with Entrez.efetch(
                db=database,
                id=accession_id,
                rettype="gb",
                retmode=retmode,
                api_key=api_key,
            ) as handle:

                for record in SeqIO.parse(handle, "genbank"):
                    seq_records.append(record)
                    pbar.update(1)

    return seq_records


# def _seqio_to_dna_info(cls, entry: SeqIO):
#     regions = []
#     for feature in entry.features:
#         if feature.type == "gene":
#             regions.append(
#                 DNARegion(
#                     name=feature.qualifiers["gene"][0],
#                     spans=[
#                         Span(
#                             start=int(feature.location.start),
#                             end=int(feature.location.end),
#                         )
#                     ],
#                     cross_reference=re.findall(
#                         GENEID_PATTERN, "".join(feature.qualifiers["db_xref"])
#                     )[0],
#                     note=feature.qualifiers["note"][0],
#                     type=DNARegionType.GENE,
#                 )
#             )

#         if feature.type == "CDS":
#             regions.append(
#                 DNARegion(
#                     name=feature.qualifiers["product"][0],
#                     spans=[
#                         Span(
#                             start=int(feature.location.start),
#                             end=int(feature.location.end),
#                         )
#                     ],
#                     cross_reference=feature.qualifiers["protein_id"][0],
#                     note=feature.qualifiers["note"][0],
#                     type=DNARegionType.CODING_SEQUENCE,
#                 )
#             )

#         if feature.type == "source":
#             organism = get_organism(entry.annotations, feature)

#     return cls(
#         source_id=entry.id,
#         name=entry.description,
#         sequence=str(entry.seq),
#         regions=regions,
#         organism=organism,
#     )
