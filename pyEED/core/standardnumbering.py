import sdRDM

from typing import List, Optional
from pydantic import Field
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

    reference_id: Optional[str] = Field(
        default=None,
        description="Standard numbering of the reference sequence",
    )

    numbered_id: Optional[str] = Field(
        default=None,
        description="Standard numbering of the query sequence",
    )

    numbering: List[str] = Field(
        description="Standard numbering of the aligned sequence",
        default_factory=ListPlus,
        multiple=True,
    )
