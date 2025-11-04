from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from itertools import combinations
from multiprocessing import cpu_count
from typing import Literal

from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    ProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich.text import Text
from skbio.alignment import pair_align as skbio_pair_align
from skbio.alignment._pair import PairAlignResult
from skbio.sequence import Sequence, SubstitutionMatrix


@dataclass
class AlignmentResult:
    """Alignment statistics and metrics for a pairwise sequence alignment.

    Attributes:
        query_id: Identifier for the query sequence.
        target_id: Identifier for the target sequence.
        score: Alignment score from the substitution matrix.
        identity: Fraction of matching positions among aligned positions.
        matches: Number of matching aligned positions.
        mismatches: Number of mismatched aligned positions.
        gaps: Number of positions with gaps in either sequence.
        min_cov: Minimum coverage between the two sequences.
        max_cov: Maximum coverage between the two sequences.
        mean_cov: Mean coverage between the two sequences.
    """

    query_id: str
    target_id: str
    score: float
    identity: float
    matches: int
    mismatches: int
    gaps: int
    min_cov: float
    max_cov: float
    mean_cov: float


def alignment_counts(seq1: str, seq2: str) -> dict[str, int]:
    """Count matches, mismatches, and gaps in aligned sequences.

    Args:
        seq1: First aligned sequence string.
        seq2: Second aligned sequence string (must have same length as seq1).

    Returns:
        Dictionary with keys "matches", "mismatches", and "gaps".
    """
    if len(seq1) != len(seq2):
        msg = "Sequences must have equal length"
        raise ValueError(msg)
    matches = sum((a != "-") and (b != "-") and (a == b) for a, b in zip(seq1, seq2, strict=False))
    mismatches = sum((a != "-") and (b not in {"-", a}) for a, b in zip(seq1, seq2, strict=False))
    gaps = sum((a == "-") or (b == "-") for a, b in zip(seq1, seq2, strict=False))
    return {"matches": matches, "mismatches": mismatches, "gaps": gaps}


def coverage(seq1: str, seq2: str) -> dict[str, float]:
    """Calculate identity and coverage metrics for aligned sequences.

    Args:
        seq1: First aligned sequence string.
        seq2: Second aligned sequence string.

    Returns:
        Dictionary with keys "identity", "min_cov", "max_cov", and "mean_cov".
    """
    counts = alignment_counts(seq1, seq2)
    aligned = counts["matches"] + counts["mismatches"]
    identity_val = 0.0 if aligned == 0 else counts["matches"] / aligned

    cov_a = aligned / len(seq1)
    cov_b = aligned / len(seq2)
    min_cov, max_cov = sorted([cov_a, cov_b])
    mean_cov = (cov_a + cov_b) / 2

    return {"identity": identity_val, "min_cov": min_cov, "max_cov": max_cov, "mean_cov": mean_cov}


def identity(seq1: str, seq2: str) -> float:
    """Calculate sequence identity as matches divided by minimum sequence length.

    Args:
        seq1: First sequence string.
        seq2: Second sequence string.

    Returns:
        Fraction of matching positions relative to the shorter sequence length.
    """
    n_aligned = alignment_counts(seq1, seq2)["matches"]
    return n_aligned / min(len(seq1), len(seq2))


def get_alignment_results(
    seq1: Sequence, seq2: Sequence, alignment: PairAlignResult
) -> AlignmentResult:
    """Convert a PairAlignResult to an AlignmentResult with computed metrics.

    Args:
        seq1: Sequence object.
        seq2: Sequence object.
        alignment: Pairwise alignment result from skbio.

    Returns:
        AlignmentResult with all computed alignment metrics.
    """
    aligned_seq1, aligned_seq2 = alignment.paths[0].to_aligned((seq1, seq2))
    counts = alignment_counts(aligned_seq1, aligned_seq2)
    cov_metrics = coverage(aligned_seq1, aligned_seq2)

    return AlignmentResult(
        query_id=seq1.metadata["id"],
        target_id=seq2.metadata["id"],
        score=alignment.score,
        identity=cov_metrics["identity"],
        matches=counts["matches"],
        mismatches=counts["mismatches"],
        gaps=counts["gaps"],
        min_cov=cov_metrics["min_cov"],
        max_cov=cov_metrics["max_cov"],
        mean_cov=cov_metrics["mean_cov"],
    )


def to_annotated_sequence(ids: list[str], sequences: list[str]) -> list[Sequence]:
    """Convert sequence IDs and sequences to annotated Sequence objects.

    Args:
        ids: List of sequence identifiers.
        sequences: List of sequence strings (must match length of ids).

    Returns:
        List of Sequence objects with metadata containing the sequence IDs.
    """
    return [
        Sequence(sequence, metadata={"id": seq_id})
        for seq_id, sequence in zip(ids, sequences, strict=False)
    ]


def build_unique_pairs(sequences: list[Sequence]) -> list[tuple[Sequence, Sequence]]:
    """Generate all unique pairwise combinations of sequences.

    Args:
        sequences: List of Sequence objects.

    Returns:
        List of tuples containing all unique pairs of sequences.
    """
    return list(combinations(sequences, 2))


def _align_task(
    seq1: Sequence,
    seq2: Sequence,
    mode: Literal["global", "local"],
    substitution_matrix: str,
    gap_cost: int,
) -> PairAlignResult:
    alignment = skbio_pair_align(seq1, seq2, mode, substitution_matrix, gap_cost)
    return get_alignment_results(seq1, seq2, alignment)


def align_pairs(
    ids: list[str],
    sequences: list[str],
    mode: Literal["global", "local"] = "global",
    n_jobs: int = -1,
    substitution_matrix: str = "blosum62",
    gap_cost: int = 2,
    verbose: bool = True,
    in_notebook: bool = False,
) -> list[PairAlignResult]:
    """Perform pairwise sequence alignment in parallel.

    Args:
        pairs: List of sequence pairs to align.
        mode: Alignment mode, either "global" or "local".
        n_jobs: Number of parallel jobs. -1 uses all available CPU cores.
        substitution_matrix: Name of substitution matrix (e.g., "blosum62").
        gap_cost: Cost for gap insertion/extension.

    Returns:
        List of PairAlignResult objects, one per input pair.
    """
    sequences = to_annotated_sequence(ids, sequences)
    pairs = build_unique_pairs(sequences)

    if n_jobs <= -1:
        n_jobs = cpu_count()

    sub_matrix = SubstitutionMatrix.by_name(substitution_matrix)
    results: list[PairAlignResult] = []

    console = Console(quiet=not verbose, force_jupyter=in_notebook)
    with Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        RateColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
        refresh_per_second=2,
    ) as progress:
        task = progress.add_task("⛓️ Aligning sequence pairs...", total=len(pairs))

        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            futures = {
                ex.submit(_align_task, seq1, seq2, mode, sub_matrix, gap_cost)
                for seq1, seq2 in pairs
            }
            for fut in as_completed(futures):
                results.append(fut.result())
                progress.update(task, advance=1)

    return results


class RateColumn(ProgressColumn):
    def render(self, task) -> Text:
        spd = task.speed  # completed/second (EMA-smoothed)
        return Text(f"{int(spd)} it/s" if spd else "- it/s")


if __name__ == "__main__":
    path = "/home/mha/projects/proteingraph/downloads/1000seq.fasta"

    # read fasta and get sequences and accessions (sequences are multiline)
    sequences = []
    accessions = []
    sequence = ""
    accession = ""
    with open(path) as f:
        for idx, line in enumerate(f):
            if line.startswith(">"):
                if idx != 0:
                    sequences.append(sequence)
                    accessions.append(accession)
                accession = line.strip().split()[0][1:]
                sequence = ""
            else:
                sequence += line.strip()
        # Don't forget last sequence
        if sequence:
            sequences.append(sequence)
            accessions.append(accession)

    print(f"Number of sequences: {len(sequences)}")
    print(f"Number of accessions: {len(accessions)}")

    results = align_pairs(accessions, sequences, n_jobs=1)
    print(f"Alignment results: {(results)}")
