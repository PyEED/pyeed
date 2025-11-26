from typing import Annotated
from uuid import uuid4

from pydantic import Field

from .annotationtype import AnnotationType
from .pyeedbase import LabelProperty, BaseNode


class Annotation(BaseNode):
    """Sequence annotation with positions and metadata."""

    id: Annotated[str, LabelProperty(unique=True, index=True)] = Field(
        default_factory=lambda: str(uuid4()),
        description="Annotation identifier",
    )
    annotation_type: AnnotationType = Field(
        ...,
        description="Type of annotation",
    )
    positions: list[int] = Field(
        ...,
        description="Sorted list of positions",
    )
    description: str | None = Field(
        default=None,
        description="Description of the annotation",
    )
