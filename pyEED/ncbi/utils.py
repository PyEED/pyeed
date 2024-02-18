from tqdm import tqdm
from time import time
from typing import List

from pyeed.core.proteininfo import ProteinInfo
from pyeed.ncbi.seq_io import get_ncbi_entrys


def load_accessions(
    accession_ids: List[str],
    database: str = "protein",
    email: str = None,
    api_key: str = None,
) -> List[ProteinInfo]:
    if not isinstance(accession_ids, list):
        try:
            import toml

            # Load the TOML file
            with open(accession_ids, "r") as toml_file:
                seed_data = toml.load(toml_file)

            # Access the list
            for value in seed_data.values():
                if isinstance(value, list):
                    accession_ids = value
        except:
            raise ValueError("Accessions must be a list or a TOML file")

    print(api_key)

    seq_entries = get_ncbi_entrys(
        accession_ids=accession_ids,
        database="protein",
        email=email,
        api_key=api_key,
    )

    protein_infos = []
    for record in seq_entries:
        protein_infos.append(ProteinInfo._from_seq_record(record))

    return protein_infos
