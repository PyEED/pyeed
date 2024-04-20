
from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator
from .abstractregion import AbstractRegion
from .dnaregiontype import DNARegionType


@forge_signature
class DNARegion(AbstractRegion):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("dnaregionINDEX"),
        xml="@id",
    )

    type: Optional[DNARegionType] = Field(
        default=None,
        description="Type of the region within the nucleotide sequence",
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="2c478e9b9618bfdc095c0c8906fbe67c80a3e2d7"
    )
