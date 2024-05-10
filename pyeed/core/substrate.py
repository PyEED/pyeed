import sdRDM

from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator
from .citation import Citation


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

    citation: Optional[Citation] = Field(
        description="Citations of the substrate",
        default_factory=Citation,
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="3607c4e340ae59061cd0b3fe9e724e58e70e0885"
    )
