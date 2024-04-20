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
        default="2c478e9b9618bfdc095c0c8906fbe67c80a3e2d7"
    )
