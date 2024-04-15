from pydantic import Field, validator
from typing import TYPE_CHECKING

from pyeed.align import AbstractAligner
from Bio.Align import PairwiseAligner as BioPairwiseAligner

if TYPE_CHECKING:
    from Bio.Align import Alignment as BioAlignment
    from Bio.Align.substitution_matrices import Array as BioSubstitutionMatrix


class PairwiseAligner(AbstractAligner):

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

    def align(self) -> "BioAlignment":
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

        alignment_result = aligner.align(self.sequences[0], self.sequences[1])[0]

        return alignment_result

    def _load_substitution_matrix(self) -> "BioSubstitutionMatrix":
        from Bio.Align import substitution_matrices

        return substitution_matrices.load(self.substitution_matrix)
