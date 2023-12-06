import sdRDM

from typing import Optional, Union, List
from pydantic import PrivateAttr, Field, validator
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .proteinregiontype import ProteinRegionType
from .abstractregion import AbstractRegion
from .span import Span


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
