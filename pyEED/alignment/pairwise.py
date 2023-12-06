from ctypes import alignment
from typing import List
from itertools import combinations
from Bio.Align import PairwiseAligner
from tqdm import tqdm

from pyEED.core import ProteinInfo
from pyEED.core import Alignment

from joblib import Parallel, delayed, cpu_count


def multi_pairwise_alignment(
    protien_infos: List[ProteinInfo],
    mode: str = "global",
    match: int = 1,
    mismatch: int = -1,
    gap_open: int = -1,
    gap_extend: int = 0,
    substitution_matrix: str = "None",
    n_jobs: int = None,
):
    pairs = list(combinations(protien_infos, 2))

    if n_jobs is None:
        n_jobs = cpu_count()

    alignments = Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(pairwise_alignment)(
            reference,
            query,
            mode,
            match,
            mismatch,
            gap_open,
            gap_extend,
            substitution_matrix,
        )
        for reference, query in tqdm(pairs, desc="⛓️ Aligning sequences")
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
) -> Alignment:
    """Creates a pairwise alignment between two sequences.

    Args:
        reference_info (ProteinInfo): Sequence 1 to be aligned.
        query_info (ProteinInfo): Sequence 2 to be aligned.
        mode (str, optional): If the alignment should be global or local. Defaults to "global".
        match (int, optional): Score of a match. Defaults to 1.
        mismatch (int, optional): Score of a mismatch. Defaults to -1.
        gap_open (int, optional): Gap open cost. Defaults to -1.
        gap_extend (int, optional): Gap extend cost. Defaults to 0.
        substitution_matrix (str, optional): Substitution matrix name. Defaults to "None".

    Raises:
        ValueError: If the mode is not global or local.
        FileNotFoundError: If the substitution matrix is not valid.

    Returns:
        PairwiseAlignment: Pairwise alignment object with the alignment scores.
    """
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
    # print(type(alignment_result))
    # print(dir(alignment_result))
    # print(alignment_result)
    print(alignment_result.format())
    gaps = alignment_result.counts().gaps
    mismatches = alignment_result.counts().mismatches
    identities = alignment_result.counts().identities

    identity = identities / len(reference_info.sequence)

    alignment = Alignment(
        reference_seq=reference_info,
        query_seqs=[query_info],
        score=alignment_result.score,
        identity=identity,
        gaps=gaps,
        mismatches=mismatches,
    )

    return alignment
