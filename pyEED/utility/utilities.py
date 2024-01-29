from typing import List

from pyEED.core.abstractsequence import AbstractSequence
from pyEED.core.sequence import Sequence


def create_multifaster(sequences: List[Sequence]) -> str:
    """Creates a multifasta string from a list of sequences."""
    return "\n".join([f">{seq.source_id}\n{seq.sequence}" for seq in sequences])
