import numpy as np

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


    def test_align_pairwise_pair_results(self):
        aligner = PairwiseAligner(mode="global")
        seq1 = "ATCG"
        seq2 = "AGTC"
        alignment_result = aligner.align_pairwise(seq1, seq2)
        expected_keys = ["seq1", "seq2", "score", "mismatches", "gaps", "identity", "start", "end"]
        assert isinstance(alignment_result, dict)
        assert all(key in alignment_result for key in expected_keys)
        assert alignment_result['seq1'] == "A-TCG"
        assert alignment_result['seq2'] == "AGTC-"
        assert alignment_result['score'] == 1.0
        assert alignment_result['mismatches'] == 1.0
        assert alignment_result['gaps'] == 0.3333333333333333
        assert alignment_result['identity'] == 0.6
        assert (alignment_result['start'] == np.array([[0, 1], [1, 3]])).all()
        assert (alignment_result['end'] == np.array([[0, 1], [2, 4]])).all()


    def test_align_multi_pairwise_results(self):
        """
        [{'seq1': 'A-TCG', 'seq2': 'AGTC-', 'score': 1.0, 'mismatches': 1.0, 'gaps': 0.3333333333333333, 'identity': 0.6, 'start': array([[0, 1],
       [1, 3]]), 'end': array([[0, 1],
       [2, 4]])}, {'seq1': '---ATCG', 'seq2': 'AAAA---', 'score': -1.0, 'mismatches': 1.0, 'gaps': 0.14285714285714285, 'identity': 0.14285714285714285, 'start': array([[0, 1]]), 'end': array([[3, 4]])}, {'seq1': 'ATCG', 'seq2': 'TTTT', 'score': -2.0, 'mismatches': 0.25, 'gaps': 1.0, 'identity': 0.25, 'start': array([[0, 4]]), 'end': array([[0, 4]])}, {'seq1': '---AGTC', 'seq2': 'AAAA---', 'score': -1.0, 'mismatches': 1.0, 'gaps': 0.14285714285714285, 'identity': 0.14285714285714285, 'start': array([[0, 1]]), 'end': array([[3, 4]])}, {'seq1': 'AGTC', 'seq2': 'TTTT', 'score': -2.0, 'mismatches': 0.25, 'gaps': 1.0, 'identity': 0.25, 'start': array([[0, 4]]), 'end': array([[0, 4]])}, {'seq1': '----AAAA', 'seq2': 'TTTT----', 'score': -2.0, 'mismatches': 1.0, 'gaps': 0.1111111111111111, 'identity': 0.0, 'start': array([], shape=(0, 2), dtype=int64), 'end': array([], shape=(0, 2), dtype=int64)}]
        """
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
        assert alignments[0]['seq1'] == "A-TCG"
        assert alignments[0]['seq2'] == "AGTC-"
        assert alignments[0]['score'] == 1.0
        assert alignments[0]['mismatches'] == 1.0
        assert alignments[0]['gaps'] == 0.3333333333333333
        assert alignments[0]['identity'] == 0.6
        assert (alignments[0]['start'] == np.array([[0, 1], [1, 3]])).all()
        assert (alignments[0]['end'] == np.array([[0, 1], [2, 4]])).all()

        assert alignments[1]['seq1'] == "---ATCG"
        assert alignments[1]['seq2'] == "AAAA---"
        assert alignments[1]['score'] == -1.0
        assert alignments[1]['mismatches'] == 1.0
        assert alignments[1]['gaps'] == 0.14285714285714285
        assert alignments[1]['identity'] == 0.14285714285714285
        assert (alignments[1]['start'] == np.array([[0, 1]])).all()
        assert (alignments[1]['end'] == np.array([[3, 4]])).all()

        assert alignments[2]['seq1'] == "ATCG"
        assert alignments[2]['seq2'] == "TTTT"
        assert alignments[2]['score'] == -2.0
        assert alignments[2]['mismatches'] == 0.25
        assert alignments[2]['gaps'] == 1.0    