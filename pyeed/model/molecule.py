from typing import Annotated, ClassVar, List, Optional

from pydantic import Field

from .embedding import Embedding
from .pyeedbase import LabelProperty, ParentReference, PyeedBase


class Molecule(PyeedBase):
    """Chemical molecule information."""

    chebi_id: Annotated[str, LabelProperty(unique=True)] = Field(
        description="ChEBI identifier",
    )
    name: Optional[str] = Field(
        None,
        description="Molecule name",
    )
    smiles: Optional[str] = Field(
        None,
        description="SMILES representation",
    )
    inchi: Optional[str] = Field(
        None,
        description="InChI representation",
    )
    embedding: List[Embedding] = Field(
        default_factory=list,
        description="Embedding vector",
    )
    PARENT_REF: ClassVar[ParentReference] = ParentReference(
        parent_node_name="Reaction",
        rel_name="PART_OF",
    )
