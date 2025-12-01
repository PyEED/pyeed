from typing import Annotated

from pydantic import Field

from .pyeedbase import BaseNode, LabelProperty


class Taxon(BaseNode):
    """Generic taxon node used for lineage and parent."""

    id: Annotated[str, LabelProperty(index=True)] = Field(
        ...,
        description="Taxonomy ID",
    )
    scientific_name: str | None = Field(
        default=None,
        description="Scientific name of the taxon",
    )
    common_name: str | None = Field(
        default=None,
        description="Common name of the taxon",
    )
    rank: str | None = Field(
        default=None,
        description="Rank of the taxon",
    )
    hidden: bool | None = Field(
        default=None,
        description="Whether the taxon is hidden",
    )
    synonyms: list[str] = Field(
        default_factory=list,
        description="Synonyms of the taxon",
    )
