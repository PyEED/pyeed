import sdRDM

from typing import Optional
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class DNASequence(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("dnasequenceINDEX"),
        xml="@id",
    )

    protein_sequence_id: str = Field(
        ...,
        description=(
            "Reference to the corresponding protein sequence to which this DNA sequence"
            " translates"
        ),
    )
