import sdRDM

from typing import Optional, Union, List, Dict
from pydantic import PrivateAttr, model_validator, validator
from uuid import uuid4
from pydantic_xml import attr, element, wrapped
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from pyeed.aligners.pairwise import PairwiseAligner
from pyeed.containers.abstract_container import AbstractContainer
from .standardnumbering import StandardNumbering
from .sequence import Sequence
from .abstractsequence import AbstractSequence


@forge_signature
class Alignment(
    sdRDM.DataModel,
    nsmap={
        "": "https://github.com/PyEED/pyeed@3b002efa6bf51e951767d8a7749ebad563897cb8#Alignment"
    },
):
    """"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    input_sequences: List[Sequence] = wrapped(
        "input_sequences",
        element(
            description="Sequences of the alignment",
            default_factory=ListPlus,
            tag="Sequence",
            json_schema_extra=dict(multiple=True),
        ),
    )

    method: Optional[str] = element(
        description="Applied alignment method",
        default=None,
        tag="method",
        json_schema_extra=dict(),
    )

    consensus: Optional[str] = element(
        description="Consensus sequence of the alignment",
        default=None,
        tag="consensus",
        json_schema_extra=dict(),
    )

    aligned_sequences: List[Sequence] = wrapped(
        "aligned_sequences",
        element(
            description="Aligned sequences of the alignment",
            default_factory=ListPlus,
            tag="Sequence",
            json_schema_extra=dict(multiple=True),
        ),
    )

    standard_numberings: List[StandardNumbering] = wrapped(
        "standard_numberings",
        element(
            description="Standard numbering of the aligned sequences",
            default_factory=ListPlus,
            tag="StandardNumbering",
            json_schema_extra=dict(multiple=True),
        ),
    )
    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="3b002efa6bf51e951767d8a7749ebad563897cb8"
    )
    _raw_xml_data: Dict = PrivateAttr(default_factory=dict)

    @model_validator(mode="after")
    def _parse_raw_xml_data(self):
        for attr, value in self:
            if isinstance(value, (ListPlus, list)) and all(
                (isinstance(i, _Element) for i in value)
            ):
                self._raw_xml_data[attr] = [elem2dict(i) for i in value]
            elif isinstance(value, _Element):
                self._raw_xml_data[attr] = elem2dict(value)
        return self

    def add_to_input_sequences(
        self,
        source_id: Optional[str] = None,
        sequence: Optional[str] = None,
        id: Optional[str] = None,
    ) -> Sequence:
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
    ) -> Sequence:
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
    ) -> StandardNumbering:
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

    @validator("input_sequences", pre=True)
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
                "Invalid sequence type. Sequences must be of type AbstractSequence or"
                " Sequence"
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

    def align(self, aligner: Union[AbstractContainer, PairwiseAligner], **kwargs):

        if issubclass(aligner, AbstractContainer):
            self._container_align(aligner, **kwargs)

        elif issubclass(aligner, PairwiseAligner):
            self._python_align(aligner, **kwargs)

        else:
            raise ValueError(
                "aligner must be an instance of AbstractContainer or PairwiseAligner"
            )

    def _container_align(self, aligner: AbstractContainer, **kwargs):
        sequences = [seq.fasta_string() for seq in self.input_sequences]

        alignment = aligner().align(sequences=sequences)

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

        self.apply_standard_numbering()

    def _python_align(self, aligner: PairwiseAligner, **kwargs):
        if len(self.input_sequences) == 2:
            alignment_reuslt = aligner(
                sequences=[
                    self.input_sequences[0].sequence,
                    self.input_sequences[1].sequence,
                ],
                **kwargs,
            ).align()

        # self.aligned_sequences = [
        #     Sequence(source_id=seq.id, sequence=str(seq.seq))
        #     for seq in alignment_reuslt
        # ]

        # TODO: Alignment has no ID attriburte
        # TODO: Map to data model
        return alignment_reuslt

    @classmethod
    def from_sequences(
        cls,
        sequences: List[AbstractSequence],
        aligner: Union[AbstractContainer, PairwiseAligner] = None,
        **kwargs,
    ):
        alignment = cls(
            input_sequences=sequences,
        )

        if aligner is not None:
            alignment.align(aligner, **kwargs)

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
