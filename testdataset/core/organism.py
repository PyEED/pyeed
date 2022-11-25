import sdRDM

from typing import Optional, Union
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Organism(sdRDM.DataModel):
    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("organismINDEX"),
        xml="@id",
    )

    ncbi_taxonomy_id: str = Field(
        ..., description="NCBI Taxonomy ID to identify the organism"
    )

    __repo__: Optional[str] = PrivateAttr(default="git://github.com/maxim945/test5.git")

    __commit__: Optional[str] = PrivateAttr(
        default="6d73cfc10f69971242fd92fe6867211fe1e1595a"
    )
