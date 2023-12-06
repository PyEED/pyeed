import sdRDM

from typing import Optional, Union, List
from pydantic import PrivateAttr, Field, validator
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class StandardNumbering(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("standardnumberingINDEX"),
        xml="@id",
    )

    sequence_id: Optional[str] = Field(
        default=None,
        description="Identifier of the aligned sequence",
    )

    numbering: List[str] = Field(
        description="Standard numbering of the aligned sequence",
        default_factory=ListPlus,
        multiple=True,
    )
