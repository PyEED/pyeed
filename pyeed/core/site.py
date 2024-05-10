import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
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

    cross_ref: Optional[str] = Field(
        default=None,
        description="Database cross reference",
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="3607c4e340ae59061cd0b3fe9e724e58e70e0885"
    )
