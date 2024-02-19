from pydantic import BaseModel, Field, validator
import sdRDM
from typing import Any, List, Optional
from itertools import combinations
from Bio.Align import PairwiseAligner as BioPairwiseAligner
from tqdm import tqdm
from pyeed.core.abstractsequence import AbstractSequence

from pyeed.core.pairwisealignment import PairwiseAlignment
from pyeed.core.sequence import Sequence
from pyeed.core import Alignment
from pyeed.core import StandardNumbering

from joblib import Parallel, delayed, cpu_count


class PairwiseAligner(BaseModel):

    sequences: List[Sequence] = Field(
        description="Sequences to be aligned",
        max_items=2,
        min_items=2,
    )
    mode: str = Field(
        description="Alignment mode",
        default="global",
    )
    match: int = Field(
        description="Score of a match",
        default=1,
    )
    mismatch: int = Field(
        description="Score of a mismatch",
        default=-1,
    )
    gap_open: int = Field(
        description="Gap open cost",
        default=-1,
    )
    gap_extend: int = Field(
        description="Gap extend cost",
        default=0,
    )
    substitution_matrix: str = Field(
        description="Substitution matrix name",
        default="None",
    )

    @validator("sequences", pre=True)
    def sequences_validator(cls, sequences):
        if all(isinstance(seq, AbstractSequence) for seq in sequences):
            return [
                Sequence(source_id=seq.source_id, sequence=seq.sequence)
                for seq in sequences
            ]
        elif all(isinstance(seq, Sequence) for seq in sequences):
            return sequences
        else:
            raise ValueError(
                "Invalid sequence type. Sequences must be of type AbstractSequence or Sequence"
            )

    @validator("mode")
    def mode_validator(cls, mode):
        modes = ["global", "local"]

        if mode not in modes:
            raise ValueError(
                f"Invalid alignment mode: {mode}. Valid modes are: {modes}"
            )

        return mode

    @validator("substitution_matrix")
    def substitution_matrix_validator(cls, substitution_matrix):
        matrices = ["BLOSUM62", "BLOSUM45", "BLOSUM80", "PAM250", "PAM30", "PAM70"]

        if substitution_matrix not in matrices:
            raise ValueError(
                f"Invalid substitution matrix: {substitution_matrix}. Available matrices are: {matrices}"
            )

        return substitution_matrix

    def align(self):
        """
        Aligns two sequences using the specified alignment parameters of the `PairwiseAligner` class.

        Returns:
            PairwiseAlignment: The aligned sequences along with alignment statistics.

        Raises:
            ValueError: If the sequences are not of type AbstractSequence or Sequence.
            ValueError: If the alignment mode is invalid.
            ValueError: If the substitution matrix name is invalid.
        """

        aligner = BioPairwiseAligner()
        aligner.mode = self.mode
        aligner.match_score = self.match
        aligner.mismatch_score = self.mismatch
        aligner.open_gap_score = self.gap_open
        aligner.extend_gap_score = self.gap_extend

        if self.substitution_matrix != "None":
            aligner.substitution_matrix = self._load_substitution_matrix()

        shorter_seq, longer_seq = sorted(self.sequences, key=lambda x: len(x.sequence))

        alignment_result = aligner.align(shorter_seq.sequence, longer_seq.sequence)[0]

        aligned_sequences = [
            Sequence(source_id=shorter_seq.source_id, sequence=alignment_result[0]),
            Sequence(source_id=longer_seq.source_id, sequence=alignment_result[1]),
        ]

        gaps = alignment_result.counts().gaps
        mismatches = alignment_result.counts().mismatches
        identities = alignment_result.counts().identities
        identity = identities / len(shorter_seq.sequence)

        standard_numbering = StandardNumbering(
            reference_id=shorter_seq.source_id,
            numbered_id=longer_seq.source_id,
            numbering=Alignment._get_numbering_string(
                shorter_seq.sequence, longer_seq.sequence
            ),
        )

        alignment = PairwiseAlignment(
            input_sequences=[shorter_seq, longer_seq],
            method=self.mode,
            aligned_sequences=aligned_sequences,
            standard_numberings=[standard_numbering],
            score=alignment_result.score,
            identity=identity,
            gaps=gaps,
            mismatches=mismatches,
        )

        return alignment

    def _load_substitution_matrix(self) -> Any:
        from Bio.Align import substitution_matrices

        return substitution_matrices.load(self.substitution_matrix)


# def multi_pairwise_alignment(
#     protien_infos: List[ProteinInfo],
#     mode: str = "global",
#     match: int = 1,
#     mismatch: int = -1,
#     gap_open: int = -1,
#     gap_extend: int = 0,
#     substitution_matrix: str = "None",
#     n_jobs: int = None,
# ):
#     pairs = list(combinations(protien_infos, 2))

#     if n_jobs is None:
#         n_jobs = cpu_count()

#     alignments = Parallel(n_jobs=n_jobs, prefer="processes")(
#         delayed(pairwise_alignment)(
#             reference,
#             query,
#             mode,
#             match,
#             mismatch,
#             gap_open,
#             gap_extend,
#             substitution_matrix,
#         )
#         for reference, query in tqdm(pairs, desc="⛓️ Aligning sequences")
#     )

#     return alignments
