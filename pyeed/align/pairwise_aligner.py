from itertools import combinations
from typing import Dict, List

from Bio.Align import Alignment as Alignment
from Bio.Align import PairwiseAligner as BioPairwiseAligner
from Bio.Align.substitution_matrices import Array as BioSubstitutionMatrix
from joblib import Parallel, cpu_count, delayed
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from pyeed.core.pairwisealignmentresult import PairwiseAlignmentResult


class PairwiseAligner:
    def __init__(self, mode: str) -> None:
        self.mode = mode
        self.match = 1
        self.mismatch = -1
        self.gap_open = -1
        self.gap_extend = 0
        self.substitution_matrix = "None"

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
    ) -> PairwiseAlignmentResult:
        """Aligns two sequences and returns the alignment results.

        Args:
            seq1 (Dict[str, str]): Sequence 1 to align. Key is the sequence ID.
            seq2 (Dict[str, str]): Sequence 2 to align. Key is the sequence ID.

        Returns:
            PairwiseAlignmentResult: The alignment results.
        """

        alignment_result = self._align(seq1, seq2)

        return self._map_alignment_results(alignment_result, seq1, seq2)

    def align_multipairwise(self, sequences: Dict[str, str], **kwargs) -> List[dict]:
        """ """
        # the sequences are stored in a dictionary, key = source.id value = sequence
        # this creates a list of all possible pairs of sequences, on this list the alignment has to happen
        pairs = list(combinations(sequences.keys(), 2))

        progress_bar = Progress(
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            BarColumn(),
            MofNCompleteColumn(),
            TextColumn("•"),
            TimeElapsedColumn(),
            TextColumn("•"),
            TimeRemainingColumn(),
        )

        alignments = []
        with progress_bar as progress:
            alignments = Parallel(n_jobs=cpu_count(), prefer="processes")(
                delayed(self.align_pairwise)(
                    {pair[0]: sequences[pair[0]]}, {pair[1]: sequences[pair[1]]}
                )
                for pair in progress.track(pairs, description="⛓️Aligning sequences...")
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
    ) -> PairwiseAlignmentResult:
        # this maps the alignment results to a dictionary
        # this dictionary is used to create the network graph

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

        # alignment_dic = {
        #     "seq1": alignment[0],
        #     "seq2": alignment[1],
        #     "seq1_id": seq1_id,
        #     "seq2_id": seq2_id,
        #     "score": alignment.score,
        #     "mismatches": 1 / (alignment.counts().mismatches + 1),
        #     "gaps": 1 / (alignment.counts().gaps + 1),
        #     "identity": identity,
        #     "start": alignment.aligned[0],
        #     "end": alignment.aligned[1],
        # }

        return result_dict

    def _load_substitution_matrix(self) -> "BioSubstitutionMatrix":
        from Bio.Align import substitution_matrices

        return substitution_matrices.load(self.substitution_matrix)


if __name__ == "__main__":
    aligner = PairwiseAligner(mode="global")
    seq1 = dict(s1="ATCG")
    seq2 = dict(s2="ATCC")
    seq3 = dict(s3="GGCC")
    seq4 = dict(s4="GGGC")

    all_seq = {**seq1, **seq2, **seq3, **seq4}

    print(aligner.align_pairwise(seq1=seq1, seq2=seq2))

    alignment_result = aligner.align_multipairwise(sequences=all_seq)

    print(alignment_result)

    print(alignment_result)
