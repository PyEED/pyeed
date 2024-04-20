
from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator
from .proteinregiontype import ProteinRegionType
from .abstractregion import AbstractRegion


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
        default="c7afb06dff889f46d0d56f3c0563e0698a525d5a"
    )
