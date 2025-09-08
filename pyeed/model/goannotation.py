from typing import Annotated, ClassVar, Optional, Tuple

from pydantic import Field

from .pyeedbase import Edge, LabelProperty, PyeedBase


class GOAnnotation(PyeedBase):
    """Gene Ontology annotation."""

    go_id: Annotated[str, LabelProperty(unique=True)] = Field(
        ...,
        description="Gene Ontology identifier",
    )
    term: Optional[str] = Field(
        None,
        description="GO term name",
    )
    definition: Optional[str] = Field(
        None,
        description="GO term definition",
    )
    EDGES: ClassVar[Tuple[Edge, ...]] = (
        Edge(parent_label="Protein", rel_name="HAS_GO_ANNOTATION"),
    )
