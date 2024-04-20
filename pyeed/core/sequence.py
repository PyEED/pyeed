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
        default="c7afb06dff889f46d0d56f3c0563e0698a525d5a"
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
