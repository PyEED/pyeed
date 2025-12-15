from typing import Annotated
from uuid import uuid4

from pydantic import Field

from .annotationcategroy import AnnotationCategory
from .pyeedbase import BaseNode, LabelProperty


class Annotation(BaseNode):
    """Sequence annotation with positions and metadata."""

    id: Annotated[str, LabelProperty(index=True)] = Field(
        default_factory=lambda: str(uuid4()),
        description="Annotation identifier",
    )
    category: AnnotationCategory = Field(
        ...,
        description="Category of the annotation",
    )
    name: str | None = Field(
        default=None,
        description="Type of annotation",
    )
    description: str | None = Field(
        default=None,
        description="Description of the annotation",
    )


if __name__ == "__main__":
    from rich import print

    # test the annotation class
    annotation = Annotation(
        category=AnnotationCategory.SITE,
        name="Site",
        description="Site of the annotation",
    )
    print(annotation)
