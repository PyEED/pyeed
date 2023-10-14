import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Site(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("siteINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the site",
    )

    type: Optional[str] = Field(
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

    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed.git")
    __commit__: Optional[str] = PrivateAttr(
        default="bc253c7b0c7f5a13a8c986328f4a6e67f2f36f3c"
    )
