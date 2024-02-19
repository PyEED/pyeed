import sdRDM

from typing import Optional
from pydantic import Field
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
