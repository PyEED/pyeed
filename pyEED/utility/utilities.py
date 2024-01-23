from typing import List

from pyEED.core.abstractsequence import AbstractSequence


def create_multifaster(sequences: List[AbstractSequence]) -> str:
    """Creates a multifasta string from a list of sequences."""
    return "\n".join([f">{seq.source_id}\n{seq.sequence}" for seq in sequences])
