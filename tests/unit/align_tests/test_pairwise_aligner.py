import numpy as np

from pyeed.core.proteinrecord import ProteinRecord
from pyeed.align.pairwise_aligner import PairwiseAligner


class TestPairwiseAligner():

    def test_align_default_parameters_returns_dictionary(self):
        aligner = PairwiseAligner(mode="global")
        seq1 = "ATCG"
        seq2 = "AGTC"
        alignment_result = aligner.align_pairwise(seq1, seq2, "seq1", "seq2")
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
        expected_keys = ["seq1", "seq2", "score", "mismatches", "gaps", "identity", "start", "end", "seq1_id", "seq2_id"]
        assert isinstance(alignments, list)
        assert all(isinstance(alignment, dict) for alignment in alignments)
        assert all(all(key in alignment for key in expected_keys) for alignment in alignments)


    def test_align_pairwise_pair_results(self):
        aligner = PairwiseAligner(mode="global")
        seq1 = "ATCG"
        seq2 = "AGTC"
        alignment_result = aligner.align_pairwise(seq1, seq2, "seq1", "seq2")
        expected_keys = ["seq1", "seq2", "score", "mismatches", "gaps", "identity", "start", "end", "seq1_id", "seq2_id"]
        assert isinstance(alignment_result, dict)
        assert all(key in alignment_result for key in expected_keys)
        assert alignment_result['seq1'] == "A-TCG"
        assert alignment_result['seq2'] == "AGTC-"
        assert alignment_result['seq1_id'] == "seq1"
        assert alignment_result['seq2_id'] == "seq2"
        assert alignment_result['score'] == 1.0
        assert alignment_result['mismatches'] == 1.0
        assert alignment_result['gaps'] == 0.3333333333333333
        assert alignment_result['identity'] == 0.6
        assert (alignment_result['start'] == np.array([[0, 1], [1, 3]])).all()
        assert (alignment_result['end'] == np.array([[0, 1], [2, 4]])).all()


    def test_align_multi_pairwise_results(self):
        aligner = PairwiseAligner(mode="global")
        sequences = {
            "seq1": "ATCG",
            "seq2": "AGTC",
            "seq3": "AAAA",
            "seq4": "TTTT"
        }
        alignments = aligner.align_multipairwise(sequences)
        expected_keys = ["seq1", "seq2", "score", "mismatches", "gaps", "identity", "start", "end", "seq1_id", "seq2_id"]
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


    def test_align_multi_pairwise_from_ncbi(self):
        mat_accessions = [
            "MBP1912539.1",
            "SEV92896.1",
            "MBO8174569.1",
            "WP_042680787.1",
            "NPA47376.1",
            "WP_167889085.1",
            "WP_048165429.1",
            "ACS90033.1",
        ]
        mats = ProteinRecord.get_ids(mat_accessions)

        aligner = PairwiseAligner(mode="global")

        sequences_data_align = {mat.id: mat.sequence for mat in mats}

        alignments = aligner.align_multipairwise(sequences_data_align)



        assert len(sequences_data_align) == 8