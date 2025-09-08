from typing import Annotated, ClassVar, Optional

from pydantic import Field

from .pyeedbase import EdgeMap, LabelProperty, PyeedBase


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
    edge_map: ClassVar[EdgeMap] = EdgeMap(
        rules={"Protein": "HAS_GO_ANNOTATION"},
    )
