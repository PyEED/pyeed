from typing import Annotated

from pydantic import Field, field_validator

from .pyeedbase import LabelProperty, PyeedBase


class Organism(PyeedBase):
    """Organism information."""

    id: Annotated[int, LabelProperty(unique=True, index=True)] = Field(
        ...,
        description="Organism identifier (NCBI taxonomy ID)",
    )
    name: str | None = Field(
        default=None,
        description="Organism name",
    )
    kingdom: str | None = Field(
        default=None,
        description="Organism kingdom",
    )
    phylum: str | None = Field(
        default=None,
        description="Organism phylum",
    )
    class_name: str | None = Field(
        default=None,
        description="Organism class name",
        alias="class",
    )
    order: str | None = Field(
        default=None,
        description="Organism order",
    )
    family: str | None = Field(
        default=None,
        description="Organism family",
    )
    genus: str | None = Field(
        default=None,
        description="Organism genus",
    )
    species: str | None = Field(
        default=None,
        description="Organism species",
    )

    @field_validator("tax_id")
    @classmethod
    def validate_tax_id(cls, v: int) -> int:
        if v <= 0:
            raise ValueError("Taxonomy ID must be positive")
        return v
