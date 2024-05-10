
from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator
from .abstractregion import AbstractRegion
from .proteinregiontype import ProteinRegionType


@forge_signature
class ProteinRegion(AbstractRegion):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("proteinregionINDEX"),
        xml="@id",
    )

    type: Optional[ProteinRegionType] = Field(
        default=None,
        description="Type of the region within the protein sequence",
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="3607c4e340ae59061cd0b3fe9e724e58e70e0885"
    )
