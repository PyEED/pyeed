from typing import Annotated

from neo4j import AsyncDriver
from pydantic import Field

from .molecule import Molecule
from .pyeedbase import BaseNode, LabelProperty


class Reaction(BaseNode):
    """Chemical reaction information."""

    id: Annotated[str, LabelProperty(index=True)] = Field(
        ...,
        description="Reaction identifier (RHEA ID)",
    )
    description: str | None = Field(
        None,
        description="Reaction description",
    )
    substrate_ids: list[str] = Field(
        default_factory=list,
        description="List of molecule identifiers (ChEBI IDs)",
    )
    product_ids: list[str] = Field(
        default_factory=list,
        description="List of molecule identifiers (ChEBI IDs)",
    )
    reversible: bool = Field(
        default=False,
        description="Whether the reaction is reversible",
    )

    async def relate_to_substrate(
        self, driver: AsyncDriver, substrate: Molecule | list[Molecule]
    ) -> None:
        """Create (self)-[HAS_SUBSTRATE]->(substrate) relationships.

        Args:
            driver: Neo4j async driver.
            substrate: Molecule or list of Molecules.
        """

        await self._relate(driver, "HAS_SUBSTRATE", substrate)

    async def relate_to_product(
        self, driver: AsyncDriver, product: Molecule | list[Molecule]
    ) -> None:
        """Create (self)-[HAS_PRODUCT]->(product) relationships.

        Args:
            driver: Neo4j async driver.
            product: Molecule or list of Molecules.
        """

        await self._relate(driver, "HAS_PRODUCT", product)
