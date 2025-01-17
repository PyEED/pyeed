import re
from collections import defaultdict
from enum import Enum
from typing import Dict, List

PATTERN_ALL = r".*"
PDB_PATTERN = r"\bpdb[a-zA-Z0-9]{8}\b|\b[a-zA-Z0-9]{4}\b"
UNIPROT_PATTERN = (
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
)
INTERPRO_PATTERN = r"IPR\d{6}"


class DBPattern(Enum):
    """Enum class for the different database id patterns."""

    UNIPROT = re.compile(UNIPROT_PATTERN)
    NCBI = re.compile(PATTERN_ALL)
    PDB = re.compile(PDB_PATTERN)
    INTERPRO = re.compile(INTERPRO_PATTERN)


class SortIDs:
    @staticmethod
    def sort(ids: List[str]) -> Dict[str, List[str]]:
        """Sorts the ids into lists based on the database patterns."""

        db_id_dict = defaultdict(list)
        for id in ids:
            if DBPattern.UNIPROT.value.match(id):
                db_id_dict[DBPattern.UNIPROT.name].append(
                    DBPattern.UNIPROT.value.match(id).group()
                )
            elif DBPattern.PDB.value.match(id):
                db_id_dict[DBPattern.PDB.name].append(
                    DBPattern.PDB.value.match(id).group()
                )
            elif DBPattern.INTERPRO.value.match(id):
                db_id_dict[DBPattern.INTERPRO.name].append(
                    DBPattern.INTERPRO.value.match(id).group()
                )
            else:
                db_id_dict[DBPattern.NCBI.name].append(
                    DBPattern.NCBI.value.match(id).group()
                )

        for key, value in db_id_dict.items():
            db_id_dict[key] = list(set(value))

        return db_id_dict
