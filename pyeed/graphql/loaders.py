from __future__ import annotations

from typing import Literal, TypeVar, overload

from neo4j import AsyncDriver
from pymilvus import AsyncMilvusClient
from strawberry.dataloader import DataLoader

from ..ingest.model import Annotation, GOAnnotation, Molecule, Protein, Reaction, Taxon
from ..ingest.model.pyeedbase import BaseNode
from ..queries import get_vectors
from .types import (
    AnnotationType,
    GOAnnotationType,
    MoleculeType,
    ProteinType,
    ReactionType,
    TaxonType,
)

T = TypeVar("T")


def get_protein_vector_loader(milvus_client: AsyncMilvusClient) -> DataLoader[str, list[float]]:
    async def load_fn(ids: list[str]) -> list[list[float]]:
        return await get_vectors(
            milvus_client=milvus_client,
            collection_name="pyeed",
            protein_ids=ids,
            vector_field_name="mean_pooling",
            as_list=True,
        )

    return DataLoader(load_fn)


@overload
def create_relationship_loader[S: BaseNode, E: BaseNode, T](
    source_cls: type[S],
    target_cls: type[E],
    target_type: type[T],
    neo4j_driver: AsyncDriver,
    *,
    rel_type: str | None = None,
    direction: Literal["out", "in", "both"] = "out",
    cardinality: Literal["one"],
) -> DataLoader[str, T | None]:
    """Create DataLoader for one-to-one relationships."""
    ...


@overload
def create_relationship_loader[S: BaseNode, E: BaseNode, T](
    source_cls: type[S],
    target_cls: type[E],
    target_type: type[T],
    neo4j_driver: AsyncDriver,
    *,
    rel_type: str | None = None,
    direction: Literal["out", "in", "both"] = "out",
    cardinality: Literal["many"],
) -> DataLoader[str, list[T]]:
    """Create DataLoader for one-to-many relationships."""
    ...


def create_relationship_loader[S: BaseNode, E: BaseNode, T](
    source_cls: type[S],
    target_cls: type[E],
    target_type: type[T],
    neo4j_driver: AsyncDriver,
    *,
    rel_type: str | None = None,
    direction: Literal["out", "in", "both"] = "out",
    cardinality: Literal["one", "many"] = "many",
) -> DataLoader[str, T | None] | DataLoader[str, list[T]]:
    """Generic DataLoader builder for Neo4j relationships.

    Args:
        source_cls: Source node class (e.g., Protein).
        target_cls: Target node class (e.g., Reaction).
        target_type: GraphQL type to convert to (e.g., ReactionType).
        neo4j_driver: Neo4j async driver.
        rel_type: Optional relationship type (e.g., "CATALYZES").
        direction: Relationship direction: "out", "in", or "both".
        cardinality: "one" for one-to-one, "many" for one-to-many.

    Returns:
        DataLoader that batches relationship lookups.
    """
    if cardinality == "one":

        async def load_fn_one(source_ids: list[str]) -> list[T | None]:
            responses = await source_cls.get_related(
                target_cls,
                driver=neo4j_driver,
                ids=source_ids,
                rel_type=rel_type,
                direction=direction,
            )
            # Map to input order, take first (or None)
            return [
                target_type(**responses[source_id][0].model_dump())
                if responses.get(source_id) and responses[source_id]
                else None
                for source_id in source_ids
            ]

        return DataLoader(load_fn_one)

    else:  # cardinality == "many"

        async def load_fn_many(source_ids: list[str]) -> list[list[T]]:
            responses = await source_cls.get_related(
                target_cls,
                driver=neo4j_driver,
                ids=source_ids,
                rel_type=rel_type,
                direction=direction,
            )
            # Map to input order, return entire list (or empty list)
            return [
                [target_type(**item.model_dump()) for item in responses.get(source_id, [])]
                for source_id in source_ids
            ]

        return DataLoader(load_fn_many)


# Convenience factory functions using the generic builder
def organism_of_protein(neo4j_driver: AsyncDriver) -> DataLoader[str, TaxonType | None]:
    """DataLoader for one-to-one: Protein -> Organism (Taxon)."""
    return create_relationship_loader(
        Protein,
        Taxon,
        TaxonType,
        neo4j_driver,
        rel_type="ORIGINATES_FROM",
        direction="out",
        cardinality="one",
    )


def reactions_of_protein(neo4j_driver: AsyncDriver) -> DataLoader[str, list[ReactionType]]:
    """DataLoader for one-to-many: Protein -> Reactions."""
    return create_relationship_loader(
        Protein,
        Reaction,
        ReactionType,
        neo4j_driver,
        rel_type="CATALYZES",
        direction="out",
        cardinality="many",
    )


def proteins_of_reaction(neo4j_driver: AsyncDriver) -> DataLoader[str, list[ProteinType]]:
    """DataLoader for one-to-many: Reaction -> Proteins (reverse of CATALYZES)."""
    return create_relationship_loader(
        Reaction,
        Protein,
        ProteinType,
        neo4j_driver,
        rel_type="CATALYZES",
        direction="in",
        cardinality="many",
    )


def substrates_of_reaction(
    neo4j_driver: AsyncDriver,
) -> DataLoader[str, list[MoleculeType]]:
    """DataLoader for one-to-many: Reaction -> Molecules with relationship type HAS_SUBSTRATE.

    Args:
        neo4j_driver: Neo4j async driver.
    """
    return create_relationship_loader(
        Reaction,
        Molecule,
        MoleculeType,
        neo4j_driver,
        rel_type="HAS_SUBSTRATE",
        direction="out",
        cardinality="many",
    )


def products_of_reaction(neo4j_driver: AsyncDriver) -> DataLoader[str, list[MoleculeType]]:
    """DataLoader for one-to-many: Reaction -> Molecules with relationship type HAS_PRODUCT.

    Args:
        neo4j_driver: Neo4j async driver.
    """
    return create_relationship_loader(
        Reaction,
        Molecule,
        MoleculeType,
        neo4j_driver,
        rel_type="HAS_PRODUCT",
        direction="out",
        cardinality="many",
    )


def go_annotations_of_protein(neo4j_driver: AsyncDriver) -> DataLoader[str, list[GOAnnotationType]]:
    """DataLoader for one-to-many: Protein -> GO Annotations."""
    return create_relationship_loader(
        Protein,
        GOAnnotation,
        GOAnnotationType,
        neo4j_driver,
        cardinality="many",
    )


def annotations_of_protein(neo4j_driver: AsyncDriver) -> DataLoader[str, list[AnnotationType]]:
    """DataLoader for one-to-many: Protein -> Annotations."""
    return create_relationship_loader(
        Protein,
        Annotation,
        AnnotationType,
        neo4j_driver,
        cardinality="many",
    )
