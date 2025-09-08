from typing import Annotated, ClassVar

from pydantic import Field

from .embedding import Embedding
from .pyeedbase import Edge, LabelProperty, PyeedBase


class Molecule(PyeedBase):
    """Chemical molecule information."""

    chebi_id: Annotated[str, LabelProperty(unique=True)] = Field(
        description="ChEBI identifier",
    )
    name: str | None = Field(
        None,
        description="Molecule name",
    )
    smiles: str | None = Field(
        None,
        description="SMILES representation",
    )
    inchi: str | None = Field(
        None,
        description="InChI representation",
    )
    embedding: list[Embedding] = Field(
        default_factory=list,
        description="Embedding vector",
    )
    EDGES: ClassVar[tuple[Edge, ...]] = (
        Edge(
            parent_label="Reaction",
            rel_name="HAS_SUBSTRATE",
            field_name="substrates",
        ),
        Edge(
            parent_label="Reaction",
            rel_name="HAS_PRODUCT",
            field_name="products",
        ),
    )
