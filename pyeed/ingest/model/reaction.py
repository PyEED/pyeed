from typing import Annotated

from pydantic import Field

from .pyeedbase import LabelProperty, PyeedBase


class Reaction(PyeedBase):
    """Chemical reaction information."""

    id: Annotated[str, LabelProperty(unique=True, index=True)] = Field(
        ...,
        description="Reaction identifier (RHEA ID)",
    )
    description: str | None = Field(
        None,
        description="Reaction description",
    )
    substrate_ids: list[str] = Field(
        default_factory=list,
        description="List of molecule identifiers (ChEBI IDs)",
    )
    product_ids: list[str] = Field(
        default_factory=list,
        description="List of molecule identifiers (ChEBI IDs)",
    )
    reversible: bool = Field(
        default=False,
        description="Whether the reaction is reversible",
    )
