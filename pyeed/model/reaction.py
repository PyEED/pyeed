from typing import Annotated, ClassVar

from pydantic import Field

from .molecule import Molecule
from .pyeedbase import Edge, LabelProperty, PyeedBase


class Reaction(PyeedBase):
    """Chemical reaction information."""

    rhea_id: Annotated[str, LabelProperty(unique=True)] = Field(
        ...,
        description="RHEA reaction identifier",
    )
    description: str | None = Field(
        None,
        description="Reaction description",
    )
    substrates: list[Molecule] = Field(
        default_factory=list,
        description="List of ChEBI identifiers",
    )
    products: list[Molecule] = Field(
        default_factory=list,
        description="List of ChEBI identifiers",
    )
    reversible: bool = Field(
        default=False,
        description="Whether the reaction is reversible",
    )
    EDGES: ClassVar[tuple[Edge, ...]] = (Edge(parent_label="Protein", rel_name="CATALYZES"),)
