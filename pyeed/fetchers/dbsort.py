import re
from enum import Enum
from typing import Dict, List
from collections import defaultdict

PATTERN_ALL = r".*"
PDB_PATTERN = r"\bpdb[a-zA-Z0-9]{8}\b|\b[a-zA-Z0-9]{4}\b"
UNIPROT_PATTERN = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
)


class DBPattern(Enum):
    """Enum class for the different database id patterns."""

    UNIPROT = re.compile(UNIPROT_PATTERN)
    NCBI = re.compile(PATTERN_ALL)
    PDB = re.compile(PDB_PATTERN)


class SortIDs:
    @staticmethod
    def sort(ids: List[str]) -> Dict[str, List[str]]:
        """Sorts the ids into lists based on the database patterns."""

        db_id_dict = defaultdict(list)
        for id in ids:
            if DBPattern.UNIPROT.value.match(id):
                db_id_dict[DBPattern.UNIPROT.name].append(id)
            elif DBPattern.PDB.value.match(id):
                db_id_dict[DBPattern.PDB.name].append(id)
            else:
                db_id_dict[DBPattern.NCBI.name].append(id)

        return db_id_dict
