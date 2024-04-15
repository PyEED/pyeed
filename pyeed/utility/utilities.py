from typing import List, Union

from pyeed.core.sequence import Sequence
from pyeed.core.abstractsequence import AbstractSequence


def create_multifaster(sequences: Union[List[Sequence], List[AbstractSequence]]) -> str:
    """Creates a multifasta string from a list of sequences."""
    return "\n".join([f">{seq.source_id}\n{seq.sequence}" for seq in sequences])
