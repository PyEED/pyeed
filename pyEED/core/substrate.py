import sdRDM

from typing import Optional
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Substrate(sdRDM.DataModel):
    """Promiscuous substrate of an enzyme 🧪"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("substrateINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the substrate",
    )

    inchi: Optional[str] = Field(
        default=None,
        description="InChI code of the substrate",
    )

    smiles: Optional[str] = Field(
        default=None,
        description="SMILES code of the substrate",
    )

    chebi_id: Optional[str] = Field(
        default=None,
        description="ChEBI ID of the substrate",
    )
