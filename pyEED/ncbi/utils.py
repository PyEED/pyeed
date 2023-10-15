from tqdm import tqdm
from typing import List


def get_nucleotide_sequences(protein_sequences: List["ProteinSequence"]):
    """Fetches the nucleotide sequences of the coding sequences of the given protein sequences.

    Args:
        protein_sequences (List[ProteinSequence]): List of protein sequences.

    Returns:
        None
    """
    for protein_sequence in tqdm(
        protein_sequences, desc="Fetching nucleotide sequences"
    ):
        protein_sequence.get_nucleotide_seq()
