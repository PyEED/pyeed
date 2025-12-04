from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from neo4j import AsyncDriver
from pymilvus import AsyncMilvusClient
from strawberry.dataloader import DataLoader

if TYPE_CHECKING:
    from .types import (
        AnnotationType,
        GOAnnotationType,
        MoleculeType,
        ProteinType,
        ReactionType,
        TaxonType,
    )


@dataclass
class GraphQLContext:
    """Injected context available in all resolvers."""

    neo4j_driver: AsyncDriver
    milvus_client: AsyncMilvusClient
    protein_vector_loader: DataLoader[str, list[float]]
    # Protein relationships
    organism_of_protein: DataLoader[str, TaxonType | None]
    reactions_of_protein: DataLoader[str, list[ReactionType]]
    go_annotations_of_protein: DataLoader[str, list[GOAnnotationType]]
    annotations_of_protein: DataLoader[str, list[AnnotationType]]
    # Reaction relationships
    proteins_of_reaction: DataLoader[str, list[ProteinType]]
    substrates_of_reaction: DataLoader[str, list[MoleculeType]]
    products_of_reaction: DataLoader[str, list[MoleculeType]]
