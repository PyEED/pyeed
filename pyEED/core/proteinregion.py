
from typing import Optional
from pydantic import Field
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
