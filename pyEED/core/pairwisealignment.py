from typing import Optional, List
from numpy import short
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator

from itertools import combinations
from Bio.Align import PairwiseAligner
from joblib import Parallel, delayed, cpu_count
from tqdm import tqdm

from pyEED.core.abstractsequence import AbstractSequence


from .alignment import Alignment


@forge_signature
class PairwiseAlignment(Alignment):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("pairwisealignmentINDEX"),
        xml="@id",
    )

    score: Optional[float] = Field(
        default=None,
        description="Alignment score",
    )

    identity: Optional[float] = Field(
        default=None,
        description="Ration of identical residues in the alignment",
    )

    similarity: Optional[float] = Field(
        default=None,
        description="Ration of similar residues in the alignment",
    )

    gaps: Optional[int] = Field(
        default=None,
        description="Number of gaps in the alignment",
    )

    mismatches: Optional[int] = Field(
        default=None,
        description="Number of mismatches in the alignment",
    )

    def align(
        self,
        mode: str = "global",
        match: int = 1,
        mismatch: int = -1,
        gap_open: int = -1,
        gap_extend: int = 0,
        substitution_matrix: str = "None",
    ) -> Alignment:
        if len(self.input_sequences) != 2:
            raise ValueError("Pairwise alignment requires exactly two sequences")

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

        alignment_result = aligner.align(
            self.input_sequences[0].sequence, self.input_sequences[1].sequence
        )[0]
        short_sequence = min(self.input_sequences, key=lambda x: len(x.sequence))

        self.gaps = alignment_result.counts().gaps
        self.mismatches = alignment_result.counts().mismatches

        identities = alignment_result.counts().identities
        self.identity = identities / len(short_sequence.sequence)

        return self
