from typing import Annotated, ClassVar, List, Optional

from pydantic import Field

from .molecule import Molecule
from .pyeedbase import LabelProperty, ParentReference, PyeedBase


class Reaction(PyeedBase):
    """Chemical reaction information."""

    rhea_id: Annotated[str, LabelProperty(unique=True)] = Field(
        ...,
        description="RHEA reaction identifier",
    )
    description: Optional[str] = Field(
        None,
        description="Reaction description",
    )
    substrates: List[Molecule] = Field(
        default_factory=list,
        description="List of ChEBI identifiers",
    )
    products: List[Molecule] = Field(
        default_factory=list,
        description="List of ChEBI identifiers",
    )
    reversible: bool = Field(
        default=False,
        description="Whether the reaction is reversible",
    )
    PARENT_REF: ClassVar[ParentReference] = ParentReference(
        parent_node_name="Protein",
        rel_name="CATALYZES",
    )
