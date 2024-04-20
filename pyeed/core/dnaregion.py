
from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator
from .dnaregiontype import DNARegionType
from .abstractregion import AbstractRegion


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
        default="c7afb06dff889f46d0d56f3c0563e0698a525d5a"
    )
