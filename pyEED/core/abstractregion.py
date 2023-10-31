import sdRDM

from typing import Optional
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class AbstractRegion(sdRDM.DataModel):
    """Annotation of a region within a sequence 🗺️"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("abstractregionINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the annotation",
    )

    start: int = Field(
        ...,
        description="Start position of the annotation",
    )

    end: int = Field(
        ...,
        description="End position of the annotation",
    )

    note: Optional[str] = Field(
        default=None,
        description="Information found in 'note' of an ncbi protein sequence entry",
    )

    cross_reference: Optional[str] = Field(
        default=None,
        description="Database cross reference",
    )
