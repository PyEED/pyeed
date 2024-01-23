from typing import List
import tempfile

from pyEED.core.abstractsequence import AbstractSequence


def create_multifaster(sequences: List[AbstractSequence]) -> str:
    multifasta = ""
    for seq in sequences:
        multifasta += f">{seq.source_id}\n{seq.sequence}\n"

    return multifasta
