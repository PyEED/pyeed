import sdRDM

from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Equivalence(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("equivalenceINDEX"),
        xml="@id",
    )

    reference_position: int = Field(
        ...,
        description="Equivalent position in the reference sequence",
    )

    sequence_position: int = Field(
        ...,
        description=(
            "Position that is equivalent to the reference sequence position that is"
            " also given"
        ),
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed.git")
    __commit__: Optional[str] = PrivateAttr(
        default="ebd06330df7dc0565be6a6c082743cf11e5cf272"
    )
