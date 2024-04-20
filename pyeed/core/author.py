import sdRDM

from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Author(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("authorINDEX"),
        xml="@id",
    )

    given_name: Optional[str] = Field(
        default=None,
        description="Given name of the author",
    )

    family_name: Optional[str] = Field(
        default=None,
        description="Family name of the author",
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="c7afb06dff889f46d0d56f3c0563e0698a525d5a"
    )
