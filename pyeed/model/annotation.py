from typing import Annotated, ClassVar, List, Optional, Tuple
from uuid import uuid4

from pydantic import Field, field_validator

from pyeed.model.annotationtype import AnnotationType
from pyeed.model.pyeedbase import Edge, LabelProperty, PyeedBase


class Annotation(PyeedBase):
    """Sequence annotation with positions and metadata."""

    id: Annotated[str, LabelProperty(unique=True)] = Field(
        default_factory=lambda: str(uuid4()),
        description="Annotation identifier",
    )
    annotation_type: AnnotationType = Field(
        ...,
        description="Type of annotation",
    )
    positions: List[int] = Field(
        ...,
        description="Sorted list of positions",
    )
    description: Optional[str] = Field(
        default=None,
        description="Description of the annotation",
    )
    EDGES: ClassVar[Tuple[Edge, ...]] = (
        Edge(parent_label="Protein", rel_name="HAS_ANNOTATION"),
    )

    @field_validator("positions")
    @classmethod
    def validate_positions(cls, v: List[int]) -> List[int]:
        if not v:
            raise ValueError("Positions cannot be empty")
        return sorted(set(v))
