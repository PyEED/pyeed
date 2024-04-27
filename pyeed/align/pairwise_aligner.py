from itertools import combinations
from typing import Dict, List, TYPE_CHECKING

from tqdm import tqdm
from joblib import Parallel, delayed, cpu_count
from Bio.Align import PairwiseAligner as BioPairwiseAligner


if TYPE_CHECKING:
    from Bio.Align import Alignment as Alignment
    from Bio.Align.substitution_matrices import Array as BioSubstitutionMatrix

class PairwiseAligner:

    def __init__(self, mode: str) -> None:
        self.mode = mode
        self.match = 1
        self.mismatch = -1
        self.gap_open = -1
        self.gap_extend = 0
        self.substitution_matrix = "None"

    def align_pairwise(self, seq1: str, seq2: str) -> dict:
        """
        This function aligns two sequences and returns the alignment results.

        Args:
            seq1 (str): The first sequence to be aligned.
            seq2 (str): The second sequence to be aligned.
        Returns:
            dict: The alignment results.
        """

        aligner = BioPairwiseAligner()
        aligner.mode = self.mode
        aligner.match_score = self.match
        aligner.mismatch_score = self.mismatch
        aligner.open_gap_score = self.gap_open
        aligner.extend_gap_score = self.gap_extend

        if self.substitution_matrix != "None":
            aligner.substitution_matrix = self._load_substitution_matrix()

        alignment_result = aligner.align(seq1, seq2)

        alignment_result_dic = self._map_alignment_results(alignment_result[0])

        return alignment_result_dic

    def align_multipairwise(self, sequences: Dict[str, str], **kwargs) -> List[dict]:
        """
        
        """
        # the sequences are stored in a dictionary, key = source.id value = sequence
        # this creates a list of all possible pairs of sequences, on this list the alignment has to happen
        pairs = list(combinations(sequences.keys(), 2))

        aligner = BioPairwiseAligner()
        aligner.mode = self.mode


        alignments = Parallel(n_jobs=cpu_count(), prefer="processes")(
            delayed(self.align_pairwise)(a[0], a[1])
            for a in tqdm(pairs, desc="⛓️ Running pairwise alignments")
        )

        return alignments

    def _map_alignment_results(self, alignment: "Alignment") -> List[dict]:
        # this maps the alignment results to a dictionary
        # this dictionary is used to create the network graph

        shorter_seq = min(alignment[0], alignment[1], key=lambda x: len(x))

        identities = alignment.counts().identities
        identity = identities / len(shorter_seq)

        alignment_dic = {
            "seq1": alignment[0],
            "seq2": alignment[1],
            "score": alignment.score,
            "mismatches": 1 / (alignment.counts().mismatches + 1),
            "gaps": 1 / (alignment.counts().gaps + 1),
            "identity": identity,
            "start": alignment.aligned[0],
            "end": alignment.aligned[1],
        }


        return alignment_dic

    def _load_substitution_matrix(self) -> "BioSubstitutionMatrix":
        from Bio.Align import substitution_matrices

        return substitution_matrices.load(self.substitution_matrix)
    
 