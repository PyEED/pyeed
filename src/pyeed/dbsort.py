import re
from collections import defaultdict
from enum import Enum
from typing import Dict, List

from loguru import logger

NCBI_NUCLEOTIDE_PATTERN = (
    r"^[A-NR-Z][0-9]{5}(?:\.[0-9]+)?$"  # Single-letter prefix with 5 digits, optional version
    r"|^[A-NR-Z]{2}[0-9]{6}(?:\.[0-9]+)?$"  # Two-letter prefix with 6 digits, optional version
    r"|^[A-NR-Z]{2}[0-9]{8}(?:\.[0-9]+)?$"  # Two-letter prefix with 8 digits, optional version
    r"|^[A-Z]{2}_[0-9]{6}\.[0-9]+$"  # RefSeq pattern with underscore
    r"|^[A-Z]{4}[0-9]{8,}(?:\.[0-9]+)?$"  # WGS/TSA pattern with optional version
)
NCBI_PROTEIN_PATTERN = (
    r"^(?!NM_|NR_|NG_|NC_)[A-Z]{2}_[0-9]{6}\.[0-9]+$"  # Exclude specific nucleotide prefixes
    r"|^[A-Z]{3}[0-9]{5}$"
    r"|^[A-Z]{3}[0-9]{7}$"
)
PDB_PATTERN = r"\bpdb[a-zA-Z0-9]{8}\b|\b[a-zA-Z0-9]{4}\b"
UNIPROT_PATTERN = (
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
)
INTERPRO_PATTERN = r"IPR\d{6}"


class DBPattern(Enum):
    """Enum class for the different database id patterns."""

    UNIPROT = re.compile(UNIPROT_PATTERN)
    NCBI_PROTEIN = re.compile(NCBI_PROTEIN_PATTERN)
    NCBI_NUCLEOTIDE = re.compile(NCBI_NUCLEOTIDE_PATTERN)
    PDB = re.compile(PDB_PATTERN)
    INTERPRO = re.compile(INTERPRO_PATTERN)


class SortIDs:
    @staticmethod
    def sort(ids: List[str]) -> Dict[str, List[str]]:
        """Sorts the ids into lists based on the database patterns."""
        db_id_dict = defaultdict(list)
        for id in ids:
            for pattern in DBPattern:
                if pattern_match := pattern.value.match(id):
                    db_id_dict[pattern.name].append(pattern_match.group())
                    break
            else:
                logger.warning(f"ID '{id}' did not match any known patterns.")

        return db_id_dict


if __name__ == "__main__":
    # Example usage
    ids = [
        "P12345",
        "Q8NBP7",
        "1XYZ",
        "IPR000001",
        "NCBI_ID_12345",
        "UNIQUEPDB1234",
        "P67890",
        "UNKNaaaaaOWN123",
        "M12345",  # Single-letter prefix with 5 digits
        "X56789",
        "AB123456",  # Two-letter prefix with 6 digits
        "AF456789",
        "AC12345678",  # Two-letter prefix with 8 digits
        "CH98765432",
        "NM_123456.1",  # RefSeq mRNA
        "NR_987654.2",  # RefSeq non-coding RNA
        "NG_012345.1",  # RefSeq genomic regions
        "NC_000001.11",  # RefSeq chromosome
        "AAAB01000001",  # Whole Genome Shotgun (WGS)
        "AAAA02000001",
        "GBAB01000001",  # Transcriptome Shotgun Assembly (TSA)
        "JADE03000001",
        "KY739715.1",
        "ABGJUB010000054.1",
    ]
    sorted_ids = SortIDs.sort(ids)
    for db, id_list in sorted_ids.items():
        print(f"{db}: {id_list}")
