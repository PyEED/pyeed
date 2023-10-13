import sdRDM

from typing import Optional
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator
from .organism import Organism


@forge_signature
class DNASequence(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("dnasequenceINDEX"),
        xml="@id",
    )

    sequence: str = Field(
        ...,
        description="The DNA sequence",
    )

    organism: Organism = Field(
        ...,
        description="Corresponding organism",
    )
