import sdRDM

from typing import Optional, Union
from pydantic import PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from pydantic import Field


@forge_signature
class Organism(sdRDM.DataModel):
    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("organismINDEX"),
        xml="@id",
    )
    ncbi_taxonomy_id: str = Field(
        ...,
        description="NCBI Taxonomy ID to identify the organism",
    )

    __repo__: Optional[str] = PrivateAttr(default="git://github.com/maxim945/test5.git")
    __commit__: Optional[str] = PrivateAttr(
        default="b3bc6040871f5464d6ad888e0e96d94d87192172"
    )
