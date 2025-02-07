from typing import List, TypeVar

T = TypeVar("T")


def chunks(lst: List[T], n: int) -> List[List[T]]:
    """Split a list into n-sized chunks."""
    return [lst[i : i + n] for i in range(0, len(lst), n)]


def to_fasta(seq: str) -> str:
    return f">query_sequence\n{seq}"


def dict_to_fasta(sequences: dict[str, str]) -> str:
    return "\n".join([f">{id}\n{seq}" for id, seq in sequences.items()])
