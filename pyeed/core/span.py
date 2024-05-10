import sdRDM

from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Span(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("spanINDEX"),
        xml="@id",
    )

    start: Optional[int] = Field(
        default=None,
        description="Start position of the span of a region",
    )

    end: Optional[int] = Field(
        default=None,
        description="End position of the span of a region",
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="3607c4e340ae59061cd0b3fe9e724e58e70e0885"
    )
