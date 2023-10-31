import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .proteinsitetype import ProteinSiteType


@forge_signature
class Site(sdRDM.DataModel):
    """Annotation of a site within a sequence 📍"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("siteINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the site",
    )

    type: Optional[ProteinSiteType] = Field(
        default=None,
        description="Type of the site",
    )

    positions: List[int] = Field(
        description="Positions of the site",
        default_factory=ListPlus,
        multiple=True,
    )

    cross_reference: Optional[str] = Field(
        default=None,
        description="Database cross reference",
    )
