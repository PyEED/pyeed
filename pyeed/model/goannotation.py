from typing import Annotated, ClassVar

from pydantic import Field

from .pyeedbase import Edge, LabelProperty, PyeedBase


class GOAnnotation(PyeedBase):
    """Gene Ontology annotation."""

    go_id: Annotated[str, LabelProperty(unique=True)] = Field(
        ...,
        description="Gene Ontology identifier",
    )
    term: str | None = Field(
        None,
        description="GO term name",
    )
    definition: str | None = Field(
        None,
        description="GO term definition",
    )
    EDGES: ClassVar[tuple[Edge, ...]] = (
        Edge(parent_label="Protein", rel_name="HAS_GO_ANNOTATION"),
    )
