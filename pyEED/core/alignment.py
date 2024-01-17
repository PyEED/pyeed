import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .organism import Organism
from .citation import Citation
from .standardnumbering import StandardNumbering
from .abstractsequence import AbstractSequence


@forge_signature
class Alignment(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("alignmentINDEX"),
        xml="@id",
    )

    reference_seq: Optional[AbstractSequence] = Field(
        default=None,
        description="Protein sequence used as reference",
        alias="reference",
    )

    query_seqs: List[AbstractSequence] = Field(
        description="Protein sequence used as query",
        default_factory=ListPlus,
        multiple=True,
    )

    method: Optional[str] = Field(
        default=None,
        description="Method used for the alignment",
    )

    consensus: Optional[str] = Field(
        default=None,
        description="Consensus sequence of the alignment",
    )

    score: Optional[float] = Field(
        default=None,
        description="Alignment score",
    )

    standard_numberings: List[StandardNumbering] = Field(
        description="Standard numbering of the aligned sequences",
        default_factory=ListPlus,
        multiple=True,
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

    def add_to_query_seqs(
        self,
        sequence: str,
        source_id: Optional[str] = None,
        name: Optional[str] = None,
        organism: Optional[Organism] = None,
        citation: Optional[Citation] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'AbstractSequence' to attribute query_seqs

        Args:
            id (str): Unique identifier of the 'AbstractSequence' object. Defaults to 'None'.
            sequence (): Sequence of the molecule.
            source_id (): Identifier of the sequence in the source database. Defaults to None
            name (): Name of the sequence. Defaults to None
            organism (): Corresponding organism. Defaults to None
            citation (): Publication of the sequence. Defaults to None
        """
        params = {
            "sequence": sequence,
            "source_id": source_id,
            "name": name,
            "organism": organism,
            "citation": citation,
        }
        if id is not None:
            params["id"] = id
        self.query_seqs.append(AbstractSequence(**params))
        return self.query_seqs[-1]

    def add_to_standard_numberings(
        self,
        sequence_id: Optional[AbstractSequence] = None,
        numbering: List[str] = ListPlus(),
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'StandardNumbering' to attribute standard_numberings

        Args:
            id (str): Unique identifier of the 'StandardNumbering' object. Defaults to 'None'.
            sequence_id (): Identifier of the aligned sequence. Defaults to None
            numbering (): Standard numbering of the aligned sequence. Defaults to ListPlus()
        """
        params = {"sequence_id": sequence_id, "numbering": numbering}
        if id is not None:
            params["id"] = id
        self.standard_numberings.append(StandardNumbering(**params))
        return self.standard_numberings[-1]
