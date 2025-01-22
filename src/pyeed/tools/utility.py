from typing import Generator


def chunks(lst: list, n: int) -> Generator[list, None, None]:
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def to_fasta(seq: str) -> str:
    return f">query_sequence\n{seq}"


def dict_to_fasta(sequences: dict[str, str]) -> str:
    return "\n".join([f">{id}\n{seq}" for id, seq in sequences.items()])
