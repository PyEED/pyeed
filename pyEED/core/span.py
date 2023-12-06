import sdRDM

from typing import Optional, Union, List
from pydantic import PrivateAttr, Field, validator
from sdRDM.base.listplus import ListPlus
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
