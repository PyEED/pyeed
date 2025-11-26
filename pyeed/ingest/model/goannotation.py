from typing import Annotated

from pydantic import Field

from .pyeedbase import BaseNode, LabelProperty


class GOAnnotation(BaseNode):
    """Gene Ontology annotation."""

    id: Annotated[str, LabelProperty(unique=True, index=True)] = Field(
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
