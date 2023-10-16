import sdRDM

from typing import Optional
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Organism(sdRDM.DataModel):
    """Description of an organism 🦠."""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("organismINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the organism",
    )

    taxonomy_id: str = Field(
        ...,
        description="NCBI Taxonomy ID to identify the organism",
    )
