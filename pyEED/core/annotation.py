import sdRDM

from typing import Optional, Union
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Annotation(sdRDM.DataModel):
    id: str = Field(
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

    function: str = Field(
        ...,
        description=(
            "Function that is found in the annotated amino acid or sub-sequence"
        ),
    )

    end_position: Optional[int] = Field(
        description=(
            "Optional end position if the annoation contains more than a single amino"
            " acid."
        ),
        default=None,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/PyEED/pyeed-data-model.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="9b23f1fd7a004c59bd50c5619397d5142f5754f0"
    )
