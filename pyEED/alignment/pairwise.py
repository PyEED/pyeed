from typing import List
from itertools import combinations
from Bio.Align import PairwiseAligner
from numpy import mat
from tqdm import tqdm

from pyeed.core import ProteinInfo
from pyeed.core import PairwiseAlignment


def multi_pairwise_alignment(
    protien_infos: List[ProteinInfo],
    mode: str = "global",
    match: int = 1,
    mismatch: int = -1,
    gap_open: int = -1,
    gap_extend: int = 0,
    substitution_matrix: str = "None",
):
    pairs = list(combinations(protien_infos, 2))

    alignments = []
    for pair in tqdm(pairs, desc=f"⁑ Aligning unique sequence pairs"):
        reference, query = pair
        alignments.append(
            pairwise_alignment(
                reference,
                query,
                mode=mode,
                gap_open=gap_open,
                gap_extend=gap_extend,
                match=match,
                mismatch=mismatch,
                substitution_matrix=substitution_matrix,
            )
        )

    return alignments


def pairwise_alignment(
    reference_info: ProteinInfo,
    query_info: ProteinInfo,
    mode: str = "global",
    match: int = 1,
    mismatch: int = -1,
    gap_open: int = -1,
    gap_extend: int = 0,
    substitution_matrix: str = "None",
):
    modes = ["global", "local"]

    if mode not in modes:
        raise ValueError(f"Invalid alignment mode: {mode}. Valid modes are: {modes}")

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend

    if substitution_matrix != "None":
        matrices = ["BLOSUM62", "BLOSUM45", "BLOSUM80", "PAM250", "PAM30", "PAM70"]
        try:
            from Bio.Align import substitution_matrices

            aligner.substitution_matrix = substitution_matrices.load(
                substitution_matrix
            )
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Invalid substitution matrix: {substitution_matrix}. Available matrices are: {matrices}"
            )

    alignment_result = aligner.align(reference_info.sequence, query_info.sequence)[0]
    gaps = alignment_result.counts().gaps
    mismatches = alignment_result.counts().mismatches
    identities = alignment_result.counts().identities

    identity = identities / len(reference_info.sequence)
    gaps_ratio = gaps / len(reference_info.sequence)
    mismatches_ratio = mismatches / len(reference_info.sequence)

    alignment = PairwiseAlignment(
        reference_seq=reference_info,
        query_seq=query_info,
        score=alignment_result.score,
        identity=identity,
        gaps=gaps,
        mismatches=mismatches,
    )

    return alignment
