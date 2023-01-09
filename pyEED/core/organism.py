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

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/PyEED/pyeed-data-model.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="9b23f1fd7a004c59bd50c5619397d5142f5754f0"
    )
