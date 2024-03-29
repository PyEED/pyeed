
from typing import Optional
from pydantic import Field
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
