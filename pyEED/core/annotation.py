import sdRDM

from typing import Optional
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Annotation(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("annotationINDEX"),
        xml="@id",
    )

    start_position: int = Field(
        ...,
        description=(
            "Start position of the annotation. A single start position without an end"
            " corresponds to a single amino acid"
        ),
    )

    end_position: int = Field(
        ...,
        description=(
            "Optional end position if the annoation contains more than a single amino"
            " acid."
        ),
    )

    note: Optional[str] = Field(
        default=None,
        description="Function that is found in the annotated amino acid or",
    )

    name: Optional[str] = Field(
        default=None,
        description="Additional note for the annotation",
    )

    db_xref: Optional[str] = Field(
        default=None,
        description="Database cross reference",
    )
