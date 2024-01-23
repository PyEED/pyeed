from typing import List, Literal, Optional
from pydantic import BaseModel
from abc import ABC, abstractmethod

from Bio.Align import PairwiseAligner as PA
from pyEED.containers.docker_interface import AbstractContainer

from pyEED.core import AbstractSequence
from pyEED.core import StandardNumbering
from pyEED.core import Alignment
from pyEED.containers.containers import DockerHubContainer


class AbstractPairwiseAligner(BaseModel, ABC):
    @abstractmethod
    def align(
        self, reference_sequence: AbstractSequence, query_sequence: AbstractSequence
    ):
        """Aligns two sequences and returns the alignment result.

        Args:
            reference_sequence (AbstractSequence): Reference sequence.
            query_sequence (AbstractSequence): Query sequence.
        """
        pass


class AbstractMultiSequenceAligner(BaseModel, ABC):
    @abstractmethod
    def _get_aligner(self, aligner: AbstractContainer):
        """Returns the aligner object.

        Args:
            aligner (DockerInterface): Docker interface object.
        """
        pass

    @abstractmethod
    def align(self, sequences: List[AbstractSequence]):
        """Aligns multiple sequences and returns the alignment result.

        Args:
            sequences (List[AbstractSequence]): List of sequences to be aligned.
        """
        pass


class ClustalOmegaAligner(AbstractMultiSequenceAligner):
    pass


class PairwiseAligner(AbstractPairwiseAligner):
    mode: Literal["global", "local"] = "global"
    match_score: int = 1
    mismatch_score: int = -1
    gap_open: int = -1
    gap_extend: int = 0
    substitution_matrix: Optional[str] = None

    def align(
        self, reference_sequence: AbstractSequence, query_sequence: AbstractSequence
    ):
        """
        Aligns two sequences and returns the alignment result.

        Args:
            reference_sequence (AbstractSequence): Reference sequence.
            query_sequence (AbstractSequence): Query sequence.

        Returns:
            Alignment: The alignment result, including the reference sequence, query sequences, alignment score, standard numberings, identity ratio, number of gaps, and number of mismatches.
        """
        aligner = PA()
        aligner.mode = self.mode
        aligner.match_score = self.match_score
        aligner.mismatch_score = self.mismatch_score
        aligner.open_gap_score = self.gap_open
        aligner.extend_gap_score = self.gap_extend

        if self.substitution_matrix:
            matrices = ["BLOSUM62", "BLOSUM45", "BLOSUM80", "PAM250", "PAM30", "PAM70"]
            try:
                from Bio.Align import substitution_matrices

                aligner.substitution_matrix = substitution_matrices.load(
                    self.substitution_matrix
                )
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Invalid substitution matrix: {self.substitution_matrix}. Available matrices are: {matrices}"
                )

        alignment_result = aligner.align(
            reference_sequence.sequence, query_sequence.sequence
        )[0]

        gaps = alignment_result.counts().gaps
        mismatches = alignment_result.counts().mismatches
        identities = alignment_result.counts().identities

        shortest_sequence = min(
            [reference_sequence, query_sequence], key=lambda x: len(x.sequence)
        )
        identity = identities / len(shortest_sequence.sequence)

        standard_number = standard_numbering(
            reference_info=reference_sequence, alignment_result=alignment_result
        )

        alignment = Alignment(
            reference_seq=reference_sequence,
            query_seqs=[query_sequence],
            score=alignment_result.score,
            standard_numberings=[standard_number],
            identity=identity,
            gaps=gaps,
            mismatches=mismatches,
        )

        return alignment


class ClustalOmegaAligner(AbstractMultiSequenceAligner):
    def align(self, sequences: List[AbstractSequence]):
        """Aligns multiple sequences and returns the alignment result.

        Args:
            sequences (List[AbstractSequence]): List of sequences to be aligned.
        """
        pass


def standard_numbering(
    reference_info: AbstractSequence, alignment_result  # return type?
) -> StandardNumbering:
    reference_seq_numbering = list(range(1, len(reference_info.sequence) + 1))
    query_seq_numbering = []
    counter_gap = 1
    counter_point = 1
    counter = 0
    pointer = reference_seq_numbering[0]

    for i in range(0, len(alignment_result[1][:])):
        if alignment_result[0][i] == "-":
            query_seq_numbering.append(
                str(reference_seq_numbering[counter] - 1) + "-" + str(counter_gap)
            )
            counter_gap = counter_gap + 1
            counter_point = 1

        elif alignment_result[0][i] == alignment_result[1][i]:
            query_seq_numbering.append(str(reference_seq_numbering[counter]))
            pointer = reference_seq_numbering[counter]
            counter = counter + 1
            counter_gap = 1
            counter_point = 1

        elif alignment_result[1][i] == "-":
            query_seq_numbering.append("")
            counter_gap = 1
            counter = counter + 1
            counter_point = 1

        else:
            query_seq_numbering.append(str(pointer) + "." + str(counter_point))
            counter_point = counter_point + 1
            counter_gap = 1

    standard_number = StandardNumbering(
        sequence_id=reference_info.id, numbering=query_seq_numbering
    )

    return standard_number
