from typing import Annotated, ClassVar, Optional

from pydantic import Field

from .pyeedbase import LabelProperty, ParentReference, PyeedBase


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
    PARENT_REF: ClassVar[ParentReference] = ParentReference(
        parent_node_name="Protein",
        rel_name="HAS_GO_ANNOTATION",
    )
