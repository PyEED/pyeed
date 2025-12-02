from typing import Annotated, Self

from neo4j import AsyncDriver
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

    @classmethod
    async def connect_lineage(cls, driver: AsyncDriver, lineage: list[Self]) -> None:
        """Connect the lineage of the taxon.

        Args:
            driver: Neo4j async driver.
            lineage: List of taxon nodes in the lineage.
        """
        pairs = [(lineage[i], lineage[i + 1]) for i in range(len(lineage) - 1)]
        await cls._bulk_relate(
            driver,
            "IS_A",
            pairs,
        )
