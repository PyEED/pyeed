import sdRDM

from typing import Optional, Union, List
from pydantic import PrivateAttr, Field, validator
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .span import Span
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
