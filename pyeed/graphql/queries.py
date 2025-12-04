from __future__ import annotations

import strawberry
from strawberry.types import Info

from pyeed.ingest.model import Reaction
from pyeed.ingest.model.protein import Protein
from pyeed.queries import get_similar_proteins_by_ids

from .context import GraphQLContext
from .filters import ProteinFilter, ReactionFilter
from .types import (
    ProteinType,
    ReactionType,
)


@strawberry.type
class Query:
    @strawberry.field(description="Get proteins by ID(s) and/or filter criteria")
    async def proteins(
        self,
        info: Info[GraphQLContext, None],
        ids: list[str] | None = None,
        filter: ProteinFilter | None = None,
    ) -> list[ProteinType]:
        """Get proteins by IDs and/or filter criteria.

        Examples:
            proteins(ids: ["P12345", "P0CW62"])
            proteins(filter: { seqLengthMin: 100, seqLengthMax: 200 })
            proteins(ids: ["P12345"], filter: { seqLengthMin: 150 })
        """
        driver = info.context.neo4j_driver

        # Extract filters and ranges directly
        exact_filters = {}
        ranges = None

        if filter:
            # Copy dict to avoid mutation
            filter_data = filter.__dict__.copy()

            # Extract range fields
            seq_len_min = filter_data.pop("seq_length_min", None)
            seq_len_max = filter_data.pop("seq_length_max", None)
            mol_min = filter_data.pop("mol_weight_min", None)
            mol_max = filter_data.pop("mol_weight_max", None)

            # Build exact filters (exclude None values)
            exact_filters = {k: v for k, v in filter_data.items() if v is not None}

            # Build ranges dict only if at least one range value is set
            range_dict = {}
            if seq_len_min is not None or seq_len_max is not None:
                range_dict["seq_length"] = (seq_len_min, seq_len_max)
            if mol_min is not None or mol_max is not None:
                range_dict["mol_weight"] = (mol_min, mol_max)

            ranges = range_dict if range_dict else None

        instances = await Protein.get_filtered(
            driver,
            ids=ids,
            ranges=ranges,
            **exact_filters,
        )
        return [ProteinType(**p.model_dump()) for p in instances]

    @strawberry.field(description="Get reactions by ID(s) and/or filter criteria")
    async def reactions(
        self,
        info: Info[GraphQLContext, None],
        ids: list[str] | None = None,
        filter: ReactionFilter | None = None,
    ) -> list[ReactionType]:
        driver = info.context.neo4j_driver

        exact_filters = {}
        if filter:
            exact_filters = {k: v for k, v in filter.__dict__.items() if v is not None}

        instances = await Reaction.get_filtered(
            driver,
            ids=ids,
            **exact_filters,
        )
        return [ReactionType(**i.model_dump()) for i in instances]

    @strawberry.field(description="Search for similar proteins by UniProt ID")
    async def proteinSimilaritySearch(
        self,
        info: Info[GraphQLContext, None],
        ids: str,
        limit: int,
    ) -> list[ProteinType]:
        driver = info.context.neo4j_driver
        milvus_client = info.context.milvus_client
        responses = await get_similar_proteins_by_ids(
            ids=ids,
            neo4j_driver=driver,
            milvus_client=milvus_client,
            collection_name="pyeed",
            limit=limit,
            vector_field_name="mean_pooling",
        )
        return [ProteinType(**r.model_dump()) for r in responses]
