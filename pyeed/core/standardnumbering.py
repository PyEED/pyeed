import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
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
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="c7afb06dff889f46d0d56f3c0563e0698a525d5a"
    )
