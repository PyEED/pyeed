import sdRDM

from typing import Optional
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Domain(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
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
