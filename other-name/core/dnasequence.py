import sdRDM

from typing import Optional
from pydantic import Field, PrivateAttr
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

    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed.git")
    __commit__: Optional[str] = PrivateAttr(
        default="bc253c7b0c7f5a13a8c986328f4a6e67f2f36f3c"
    )
