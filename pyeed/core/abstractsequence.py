import sdRDM

import os
from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator
from .citation import Citation
from .organism import Organism


@forge_signature
class AbstractSequence(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("abstractsequenceINDEX"),
        xml="@id",
    )

    source_id: Optional[str] = Field(
        default=None,
        description="Identifier of the sequence in the source database",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the sequence",
    )

    sequence: str = Field(
        ...,
        description="Sequence of the molecule",
    )

    organism: Optional[Organism] = Field(
        default=None,
        description="Corresponding organism",
    )

    citation: Optional[Citation] = Field(
        description="Publication of the sequence",
        default_factory=Citation,
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="2c478e9b9618bfdc095c0c8906fbe67c80a3e2d7"
    )

    def _fasta_string(self):
        return f">{self.source_id}\n{self.sequence}"

    def to_fasta(self, outdir: Optional[str] = None):
        if not outdir:
            outdir = os.getcwd() + "/" + self.source_id + ".fasta"

        with open(outdir, "w") as f:
            f.write(self._fasta_string())
