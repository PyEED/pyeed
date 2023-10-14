import sdRDM

from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Region(sdRDM.DataModel):
    """Annotation of a protein sequence."""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("regionINDEX"),
        xml="@id",
    )

    start: int = Field(
        ...,
        description=(
            "Start position of the annotation. A single start position without an end"
            " corresponds to a single amino acid"
        ),
    )

    end: int = Field(
        ...,
        description=(
            "Optional end position if the annotation contains more than a single amino"
            " acid"
        ),
    )

    note: Optional[str] = Field(
        default=None,
        description="Information found in 'note' of an ncbi protein sequence entry",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the annotation",
    )

    cross_reference: Optional[str] = Field(
        default=None,
        description="Database cross reference",
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed.git")
    __commit__: Optional[str] = PrivateAttr(
        default="73fe06fa32d40a3a1bd71455bc2052c013d24984"
    )
