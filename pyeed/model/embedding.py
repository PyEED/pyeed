import re
from typing import Annotated, Any, ClassVar, List, Tuple
from uuid import uuid4

from pydantic import Field, ValidationInfo, field_validator

from .pyeedbase import Edge, LabelProperty, PyeedBase, _is_neo4j_prop_value


class Embedding(PyeedBase):
    """Metadata about an embedding."""

    model_name: str = Field(
        ...,
        description="Name of the embedding model",
    )
    layer_index: int = Field(
        default=-1,
        description="Model layer number embedding matrix was extracted prior to pooling.",
    )
    pooling_method: str = Field(
        ...,
        description="Pooling method",
    )
    vector: Annotated[
        List[float],
        LabelProperty(vector_index=True),
    ] = Field(
        ...,
        description="The embedding vector",
    )
    n_dims: int = Field(
        ...,
        description="Embedding vector length",
    )
    id: Annotated[str, LabelProperty(unique=True)] = Field(
        default_factory=lambda: str(uuid4()),
        description="Embedding identifier",
    )
    EDGES: ClassVar[Tuple[Edge, ...]] = (
        Edge(parent_label="Protein", rel_name="HAS_EMBEDDING"),
        Edge(parent_label="Molecule", rel_name="HAS_EMBEDDING"),
    )

    @field_validator("model_name")
    @classmethod
    def model_name_slug(cls, v: str) -> str:
        # Convert to lowercase and replace non [a-z0-9_] with underscores
        v_clean = re.sub(r"[^a-z0-9_]", "_", v.lower())
        return v_clean

    @field_validator("pooling_method")
    @classmethod
    def pooling_method_slug(cls, v: str) -> str:
        # Convert to lowercase and replace non [a-z0-9_] with underscores
        v_clean = re.sub(r"[^a-z0-9_]", "_", v.lower())
        return v_clean

    @field_validator("vector")
    @classmethod
    def _check_vector_len(cls, v: List[float], info: ValidationInfo) -> List[float]:
        # Access n_dims via info.data (other fields that have already been validated)
        n_dims = info.data.get("n_dims")
        if n_dims is not None and n_dims != len(v):
            raise ValueError(f"n_dims={n_dims} != len(vector)={len(v)}")
        return v

    @property
    def neo4j_vector_prop(self) -> str:
        """Dynamic property to write into Neo4j (one ANN index per property)"""
        return f"vec__{self.model_name}__{self.pooling_method}"

    def to_dict(self) -> dict[str, Any]:
        d = self.model_dump(exclude_none=True, exclude_unset=True)
        custom = d.pop("custom", {}) or {}
        vec = d.pop("vector")
        d[self.neo4j_vector_prop] = vec
        flat = {**d, **custom}
        return {k: v for k, v in flat.items() if _is_neo4j_prop_value(v)}
