from itertools import combinations
from typing import Dict, List

from Bio.Align import Alignment as Alignment
from Bio.Align import PairwiseAligner as BioPairwiseAligner
from Bio.Align.substitution_matrices import Array as BioSubstitutionMatrix
from joblib import Parallel, cpu_count, delayed
from rich.progress import Progress


class PairwiseAligner:
    def __init__(
        self,
        mode: str = "global",
        match: int = 1,
        mismatch: int = -1,
        gap_open: int = -1,
        gap_exted: int = 0,
        substitution_matrix: str = "None",
    ) -> None:
        self.mode = mode
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extend = gap_exted
        self.substitution_matrix = substitution_matrix

    def _align(
        self,
        seq1: Dict[str, str],
        seq2: Dict[str, str],
    ) -> Alignment:
        """Aligns two sequences and returns the alignment results.

        Args:
            seq1 (Dict[str, str]): Sequence 1 to align. Key is the sequence ID.
            seq2 (Dict[str, str]): Sequence 2 to align. Key is the sequence ID.
            progress (Progress, optional): `Rich` progress. Defaults to None.
            task_id (TaskID, optional): `Rich` task_id. Defaults to None.

        Returns:
            Alignment: The alignment results.
        """

        aligner = self._get_aligner()

        results = aligner.align(list(seq1.values())[0], list(seq2.values())[0])

        return results[0]

    def align_pairwise(
        self,
        seq1: Dict[str, str],
        seq2: Dict[str, str],
    ) -> dict:
        """Aligns two sequences and returns the alignment results.

        Args:
            seq1 (Dict[str, str]): Sequence 1 to align. Key is the sequence ID.
            seq2 (Dict[str, str]): Sequence 2 to align. Key is the sequence ID.

        Returns:
            dict: Has the same signature as a `PairwiseAlignmentResult` object.
        """

        alignment_result = self._align(seq1, seq2)

        return self._map_alignment_results(alignment_result, seq1, seq2)

    def align_multipairwise(self, sequences: Dict[str, str], **kwargs) -> List[dict]:
        """Creates all possible pairwise alignments from a dictionary of sequences.

        Args:
            sequences (Dict[str, str]): A dictionary of sequences to align. The key is the sequence ID.

        Returns:
            List[dict]: A list of dictionaries containing the alignment results.
        """

        pairs = list(combinations(sequences.keys(), 2))

        with Progress() as progress:
            alignments = Parallel(n_jobs=cpu_count(), prefer="processes")(
                delayed(self.align_pairwise)(
                    {pair[0]: sequences[pair[0]]}, {pair[1]: sequences[pair[1]]}
                )
                for pair in progress.track(
                    pairs, description=f"⛓️ Aligning {len(pairs)} sequence pairs..."
                )
            )

        return alignments

    def _get_aligner(self) -> BioPairwiseAligner:
        """Creates a BioPython pairwise aligner object with the specified parameters
        from the class instance."""
        aligner = BioPairwiseAligner()
        aligner.mode = self.mode
        aligner.match_score = self.match
        aligner.mismatch_score = self.mismatch
        aligner.open_gap_score = self.gap_open
        aligner.extend_gap_score = self.gap_extend

        if self.substitution_matrix != "None":
            aligner.substitution_matrix = self._load_substitution_matrix()

        return aligner

    def _map_alignment_results(
        self, alignment: Alignment, seq1: Dict[str, str], seq2: Dict[str, str]
    ) -> dict:
        """Maps the alignment results to a dictionary.
        The dictionaly has the same signature as a `PairwiseAlignmentResult` object.

        Returns:
            dict: A dictionary containing the alignment results.
        """

        shorter_seq = min(alignment[0], alignment[1], key=lambda x: len(x))

        identities = alignment.counts().identities
        identity = identities / len(shorter_seq)
        gaps = alignment.counts().gaps
        mismatches = alignment.counts().mismatches

        sequences = [
            {"id": list(seq1.keys())[0], "sequence": list(seq1.values())[0]},
            {"id": list(seq2.keys())[0], "sequence": list(seq2.values())[0]},
        ]

        aligned_sequences = [
            {"id": list(seq1.keys())[0], "sequence": alignment[0]},
            {"id": list(seq2.keys())[0], "sequence": alignment[1]},
        ]

        result_dict = {
            "score": alignment.score,
            "identity": identity,
            "gaps": gaps,
            "mismatches": mismatches,
            "sequences": sequences,
            "aligned_sequences": aligned_sequences,
        }

        return result_dict

    def _load_substitution_matrix(self) -> "BioSubstitutionMatrix":
        from Bio.Align import substitution_matrices

        return substitution_matrices.load(self.substitution_matrix)
