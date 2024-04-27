from typing import Dict, List, TYPE_CHECKING

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

        alignment_result_dic = self._map_alignment_results(alignment_result)

        return alignment_result_dic

    def _map_alignment_results(self, alignment_results: List["Alignment"]) -> List[tuple]:
        # this maps the alignment results to a dictionary
        # this dictionary is used to create the network graph
        alignment_results_dic = []

        for alignment in alignment_results:
            shorter_seq = min(alignment[0], alignment[1], key=lambda x: len(x))

            identities = alignment.counts().identities
            identity = identities / len(shorter_seq.sequence)

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

            alignment_results_dic.append(alignment_dic)

        return alignment_results_dic

    

    def _load_substitution_matrix(self) -> "BioSubstitutionMatrix":
        from Bio.Align import substitution_matrices

        return substitution_matrices.load(self.substitution_matrix)

    def _create_pairwise_alignments(
        self, input_sequences, aligner: "PairwiseAligner", **kwargs
    ):
        """
        Creates pairwise alignments between sequences.

        This method creates pairwise alignments between sequences in the network.
        The pairwise alignments are stored in the 'pairwise_alignments' attribute of the SequenceNetwork object.
        This is done for the later visualization of the network graph with cytoscope.

        Args:
            aligner (PairwiseAligner): Python-based aligner to be called.

        Raises:
            ValueError: If the number of sequences is less than 2.

        Returns:
            Nothing the data is stored internally in fields of the class.
        """

        # Pairwise alignment
        if len(input_sequences) == 2:
            pairwise_aligner = aligner(
                sequences=[
                    input_sequences[0].sequence,
                    input_sequences[1].sequence,
                ],
                **kwargs,
            )
            alignment_result = pairwise_aligner.align()

            return self._map_pairwise_alignment_results(
                alignment_result,
                pair=(
                    input_sequences[0],
                    input_sequences[1],
                ),
                mode=pairwise_aligner.mode,
            )

        # Multi pairwise alignment
        elif len(input_sequences) > 2:
            pairs = list(combinations(input_sequences, 2))

            aligners = [
                aligner(sequences=[s.sequence for s in pair], **kwargs)
                for pair in pairs
            ]

            alignments = Parallel(n_jobs=cpu_count(), prefer="processes")(
                delayed(a.align)()
                for a in tqdm(aligners, desc="⛓️ Running pairwise alignments")
            )

            return alignments, pairs, aligners[0].mode

        else:
            raise ValueError(
                f"Alignment Error. Recieved {len(input_sequences)} sequences. Expected 2."
            )

    def align_multipairwise(self, sequences: List[Dict[str, str]]) -> List[dict]:
        pass
 