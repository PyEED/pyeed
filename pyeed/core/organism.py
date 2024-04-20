import sdRDM

from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Organism(sdRDM.DataModel):
    """Description of an organism 🦠"""

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

    domain: Optional[str] = Field(
        default=None,
        description="Domain of the organism",
    )

    kingdom: Optional[str] = Field(
        default=None,
        description="Kingdom of the organism",
    )

    phylum: Optional[str] = Field(
        default=None,
        description="Phylum of the organism",
    )

    tax_class: Optional[str] = Field(
        default=None,
        description="Class of the organism",
    )

    order: Optional[str] = Field(
        default=None,
        description="Order of the organism",
    )

    family: Optional[str] = Field(
        default=None,
        description="Family of the organism",
    )

    genus: Optional[str] = Field(
        default=None,
        description="Genus of the organism",
    )

    species: Optional[str] = Field(
        default=None,
        description="Species of the organism",
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="c7afb06dff889f46d0d56f3c0563e0698a525d5a"
    )
