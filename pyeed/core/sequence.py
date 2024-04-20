import sdRDM

from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Sequence(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("sequenceINDEX"),
        xml="@id",
    )

    source_id: Optional[str] = Field(
        default=None,
        description="Identifier of the sequence in the source database",
    )

    sequence: Optional[str] = Field(
        default=None,
        description="Sequence of the alignment. Gaps are represented by '-'",
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="2c478e9b9618bfdc095c0c8906fbe67c80a3e2d7"
    )

    def fasta_string(self) -> str:
        """
        This method returns the sequence in FASTA format

        Returns:
            str: Sequence in FASTA format
        """
        return f">{self.source_id}\n{self.sequence}"

    def __str__(self) -> str:
        return self.fasta_string()
