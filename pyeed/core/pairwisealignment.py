
from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator
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
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="c7afb06dff889f46d0d56f3c0563e0698a525d5a"
    )
