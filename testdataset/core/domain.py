import sdRDM

from typing import Optional, Union
from pydantic import PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from pydantic import Field


@forge_signature
class Domain(sdRDM.DataModel):
    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("domainINDEX"),
        xml="@id",
    )
    name: str = Field(
        ...,
        description="Name of the annotated domain",
    )

    start_position: int = Field(
        ...,
        description="Position in the sequence where the domain starts",
    )

    end_position: int = Field(
        ...,
        description="Position in the sequence where the domain ends",
    )

    __repo__: Optional[str] = PrivateAttr(default="git://github.com/maxim945/test5.git")
    __commit__: Optional[str] = PrivateAttr(
        default="3446530b04f1864b8ceaff1b5ff80f2261deecc8"
    )
