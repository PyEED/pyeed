from typing import Annotated, ClassVar, Optional

from pydantic import Field, field_validator

from .pyeedbase import EdgeMap, LabelProperty, PyeedBase


class Organism(PyeedBase):
    """Organism information."""

    tax_id: Annotated[int, LabelProperty(unique=True)] = Field(
        ...,
        description="NCBI taxonomy ID",
    )
    name: Optional[str] = Field(
        default=None,
        description="Organism name",
    )

    edge_map: ClassVar[EdgeMap] = EdgeMap(
        rules={"Protein": "ORIGINATES_FROM"},
    )

    @field_validator("tax_id")
    @classmethod
    def validate_tax_id(cls, v: int) -> int:
        if v <= 0:
            raise ValueError("Taxonomy ID must be positive")
        return v
