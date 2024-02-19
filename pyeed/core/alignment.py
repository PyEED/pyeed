import re
import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from .abstractsequence import AbstractSequence
from .sequence import Sequence
from .standardnumbering import StandardNumbering

from pyeed.containers.abstract_container import AbstractContainer


@forge_signature
class Alignment(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("alignmentINDEX"),
        xml="@id",
    )

    method: Optional[str] = Field(
        default=None,
        description="Applied alignment method",
    )

    consensus: Optional[str] = Field(
        default=None,
        description="Consensus sequence of the alignment",
    )

    input_sequences: List[Sequence] = Field(
        description="Sequences of the alignment",
        default_factory=ListPlus,
        multiple=True,
    )

    aligned_sequences: List[Sequence] = Field(
        description="Aligned sequences of the alignment",
        default_factory=ListPlus,
        multiple=True,
    )

    standard_numberings: List[StandardNumbering] = Field(
        description="Standard numbering of the aligned sequences",
        default_factory=ListPlus,
        multiple=True,
    )

    def add_to_input_sequences(
        self,
        source_id: Optional[str] = None,
        sequence: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Sequence' to attribute input_sequences

        Args:
            id (str): Unique identifier of the 'Sequence' object. Defaults to 'None'.
            source_id (): Identifier of the sequence in the source database. Defaults to None
            sequence (): Sequence of the alignment. Gaps are represented by '-'. Defaults to None
        """
        params = {"source_id": source_id, "sequence": sequence}
        if id is not None:
            params["id"] = id
        self.input_sequences.append(Sequence(**params))
        return self.input_sequences[-1]

    def add_to_aligned_sequences(
        self,
        source_id: Optional[str] = None,
        sequence: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Sequence' to attribute aligned_sequences

        Args:
            id (str): Unique identifier of the 'Sequence' object. Defaults to 'None'.
            source_id (): Identifier of the sequence in the source database. Defaults to None
            sequence (): Sequence of the alignment. Gaps are represented by '-'. Defaults to None
        """
        params = {"source_id": source_id, "sequence": sequence}
        if id is not None:
            params["id"] = id
        self.aligned_sequences.append(Sequence(**params))
        return self.aligned_sequences[-1]

    def add_to_standard_numberings(
        self,
        reference_id: Optional[str] = None,
        numbered_id: Optional[str] = None,
        numbering: List[str] = ListPlus(),
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'StandardNumbering' to attribute standard_numberings

        Args:
            id (str): Unique identifier of the 'StandardNumbering' object. Defaults to 'None'.
            reference_id (): Standard numbering of the reference sequence. Defaults to None
            numbered_id (): Standard numbering of the query sequence. Defaults to None
            numbering (): Standard numbering of the aligned sequence. Defaults to ListPlus()
        """
        params = {
            "reference_id": reference_id,
            "numbered_id": numbered_id,
            "numbering": numbering,
        }
        if id is not None:
            params["id"] = id
        self.standard_numberings.append(StandardNumbering(**params))
        return self.standard_numberings[-1]

    def align(self, aligner: AbstractContainer, **kwargs):
        alignment = aligner().align(sequences=self.input_sequences)
        self.method = aligner._container_info.name

        self.aligned_sequences = [
            Sequence(source_id=seq.id, sequence=str(seq.seq)) for seq in alignment
        ]

        # TODO: This is a workaround for the BioPython bug that causes the dumb_consensus() method to raise a warning
        # Need to find new way to calculate consensus
        import warnings
        from Bio.Align import AlignInfo

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.consensus = str(AlignInfo.SummaryInfo(alignment).dumb_consensus())

    @classmethod
    def from_sequences(cls, sequences: List[AbstractSequence]):
        alignment = cls()
        alignment.input_sequences = [
            Sequence(source_id=seq.source_id, sequence=str(seq.sequence))
            for seq in sequences
        ]
        return alignment

    @staticmethod
    def _get_numbering_string(reference: str, query: str) -> List[str]:
        """
        Assigns pairwise numbering to the reference and query sequences.

        Args:
            reference (str): The reference sequence.
            query (str): The query sequence.

        Returns:
            List[str]: A list of pairwise numbering.

        """

        numbering = []
        reference_counter = 0
        query_counter = 1

        for ref_pos, que_pos in zip(reference, query):
            if ref_pos == "-":
                numbering.append(f"{reference_counter}.{query_counter}")
                query_counter += 1
            else:
                reference_counter += 1
                if que_pos != "-":
                    numbering.append(f"{reference_counter}.{query_counter}")
                    query_counter += 1
                else:
                    numbering.append(str(reference_counter))

        return numbering

    def apply_standard_numbering(
        self,
        reference: Sequence = None,
    ):
        """
        Apply standard numbering to the aligned sequences.

        Args:
            reference (Sequence, optional): The reference sequence to use for numbering.
            If not provided, the first aligned sequence will be used as the reference.
            Defaults to None.

        Raises:
            ValueError: If the sequences are not aligned.

        """
        if not self.aligned_sequences:
            raise ValueError(
                "Sequences must be aligned first. Run the align() method first."
            )

        if reference == None:
            reference = self.aligned_sequences[0]
            aligned_sequences = self.aligned_sequences[1:]

        standard_numberings = []
        for aligned_sequence in aligned_sequences:
            numbering = self._get_numbering_string(
                reference=reference.sequence, query=aligned_sequence.sequence
            )

            standard_numberings.append(
                StandardNumbering(
                    reference_id=reference.source_id,
                    numbered_id=aligned_sequence.source_id,
                    numbering=numbering,
                )
            )

        self.standard_numberings = standard_numberings
