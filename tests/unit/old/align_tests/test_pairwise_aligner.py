# from pyeed.align.pairwise import PairwiseAligner
# from pyeed.core.proteinrecord import ProteinRecord


# class TestPairwiseAligner:
#     def test_align_default_parameters_returns_dictionary(self):
#         aligner = PairwiseAligner(mode="global")
#         seq1 = "ATCG"
#         seq2 = "AGTC"
#         alignment_result = aligner.align_pairwise(
#             seq1=dict(seq1="ATCG"), seq2=dict(seq2="AGTC")
#         )
#         expected_keys = [
#             "sequences",
#             "aligned_sequences",
#             "score",
#             "mismatches",
#             "gaps",
#             "identity",
#         ]
#         assert isinstance(alignment_result, dict)
#         assert all(key in alignment_result for key in expected_keys)

#     def test_align_multipairwise_with_more_than_two_sequences(self):
#         aligner = PairwiseAligner(mode="global")
#         sequences = {"seq1": "ATCG", "seq2": "AGTC", "seq3": "AAAA", "seq4": "TTTT"}
#         alignments = aligner.align_multipairwise(sequences)
#         assert isinstance(alignments, list)
#         assert all(isinstance(alignment, dict) for alignment in alignments)

#     def test_align_pairwise_pair_results(self):
#         aligner = PairwiseAligner(mode="global")
#         seq1 = "ATCG"
#         seq2 = "AGTC"
#         alignment_result = aligner.align_pairwise(
#             seq1=dict(seq1="ATCG"), seq2=dict(seq2="AGTC")
#         )
#         expected_keys = [
#             "sequences",
#             "aligned_sequences",
#             "score",
#             "mismatches",
#             "gaps",
#             "identity",
#         ]
#         assert isinstance(alignment_result, dict)
#         assert all(key in alignment_result for key in expected_keys)
#         assert alignment_result["aligned_sequences"][0]["sequence"] == "A-TCG"
#         assert alignment_result["aligned_sequences"][1]["sequence"] == "AGTC-"
#         assert alignment_result["score"] == 1.0
#         assert alignment_result["mismatches"] == 0.0
#         assert alignment_result["gaps"] == 2
#         assert alignment_result["identity"] == 0.6

#     def test_align_multi_pairwise_results(self):
#         aligner = PairwiseAligner(mode="global")
#         sequences = {"seq1": "ATCG", "seq2": "AGTC", "seq3": "AAAA", "seq4": "TTTT"}
#         alignments = aligner.align_multipairwise(sequences)
#         expected_keys = [
#             "sequences",
#             "aligned_sequences",
#             "score",
#             "mismatches",
#             "gaps",
#             "identity",
#         ]
#         assert isinstance(alignments, list)
#         assert all(isinstance(alignment, dict) for alignment in alignments)
#         assert all(
#             all(key in alignment for key in expected_keys) for alignment in alignments
#         )
#         assert alignments[0]["aligned_sequences"][0]["sequence"] == "A-TCG"
#         assert alignments[0]["aligned_sequences"][1]["sequence"] == "AGTC-"
#         assert alignments[0]["score"] == 1.0
#         assert alignments[0]["mismatches"] == 0
#         assert alignments[0]["gaps"] == 2
#         assert alignments[0]["identity"] == 0.6

#         assert alignments[1]["aligned_sequences"][0]["sequence"] == "---ATCG"
#         assert alignments[1]["aligned_sequences"][1]["sequence"] == "AAAA---"
#         assert alignments[1]["score"] == -1.0
#         assert alignments[1]["mismatches"] == 0
#         assert alignments[1]["gaps"] == 6
#         assert alignments[1]["identity"] == 0.14285714285714285

#         assert alignments[2]["aligned_sequences"][0]["sequence"] == "ATCG"
#         assert alignments[2]["aligned_sequences"][1]["sequence"] == "TTTT"
#         assert alignments[2]["score"] == -2.0
#         assert alignments[2]["mismatches"] == 3
#         assert alignments[2]["gaps"] == 0

#     def test_align_multi_pairwise_from_ncbi(self):
#         mat_accessions = [
#             "MBP1912539.1",
#             "SEV92896.1",
#         ]
#         mats = ProteinRecord.get_ids(mat_accessions)

#         aligner = PairwiseAligner(mode="global")

#         sequences_data_align = {mat.id: mat.sequence for mat in mats}

#         alignments = aligner.align_multipairwise(sequences_data_align)

#         expected_keys = [
#             "sequences",
#             "aligned_sequences",
#             "score",
#             "mismatches",
#             "gaps",
#             "identity",
#         ]
#         assert isinstance(alignments, list)
#         assert all(isinstance(alignment, dict) for alignment in alignments)
#         assert all(
#             all(key in alignment for key in expected_keys) for alignment in alignments
#         )
#         assert (
#             alignments[0]["aligned_sequences"][0]["id"] == "MBP1912539.1"
#             or alignments[0]["aligned_sequences"][1]["id"] == "MBP1912539.1"
#         )
#         assert (
#             alignments[0]["aligned_sequences"][0]["id"] == "SEV92896.1"
#             or alignments[0]["aligned_sequences"][1]["id"] == "SEV92896.1"
#         )
#         assert alignments[0]["score"] == 327.0
