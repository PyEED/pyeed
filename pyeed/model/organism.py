from typing import Annotated, ClassVar

from pydantic import Field, field_validator

from .pyeedbase import Edge, LabelProperty, PyeedBase


class Organism(PyeedBase):
    """Organism information."""

    tax_id: Annotated[int, LabelProperty(unique=True)] = Field(
        ...,
        description="NCBI taxonomy ID",
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

    EDGES: ClassVar[tuple[Edge, ...]] = (Edge(parent_label="Protein", rel_name="ORIGINATES_FROM"),)

    @field_validator("tax_id")
    @classmethod
    def validate_tax_id(cls, v: int) -> int:
        if v <= 0:
            raise ValueError("Taxonomy ID must be positive")
        return v
