from __future__ import annotations

import strawberry
from strawberry.types import Info

from pyeed.ingest.model.protein import Protein

from .context import GraphQLContext
from .types import (
    ProteinType,
)


@strawberry.type
class Query:
    @strawberry.field(description="Get a protein by UniProt ID")
    async def protein(
        self,
        info: Info[GraphQLContext, None],
        id: str,
    ) -> ProteinType | None:
        instance = await Protein.get(
            info.context.neo4j_driver,
            id=id,
        )
        if not instance:
            return None
        return ProteinType(**instance.model_dump())
