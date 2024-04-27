from pyeed.align.pairwise_aligner import PairwiseAligner


class TestPairwiseAligner():

    def test_align_default_parameters_returns_dictionary(self):
        aligner = PairwiseAligner(mode="global")
        seq1 = "ATCG"
        seq2 = "AGTC"
        alignment_result = aligner.align_pairwise(seq1, seq2)
        expected_keys = ["seq1", "seq2", "score", "mismatches", "gaps", "identity", "start", "end"]
        assert isinstance(alignment_result, dict)
        assert all(key in alignment_result for key in expected_keys)

    def test_align_multipairwise_with_more_than_two_sequences(self):
        aligner = PairwiseAligner(mode="global")
        sequences = {
            "seq1": "ATCG",
            "seq2": "AGTC",
            "seq3": "AAAA",
            "seq4": "TTTT"
        }
        alignments = aligner.align_multipairwise(sequences)
        expected_keys = ["seq1", "seq2", "score", "mismatches", "gaps", "identity", "start", "end"]
        assert isinstance(alignments, list)
        assert all(isinstance(alignment, dict) for alignment in alignments)
        assert all(all(key in alignment for key in expected_keys) for alignment in alignments)