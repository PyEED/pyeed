from __future__ import annotations

from typing import TYPE_CHECKING, Annotated

import strawberry
from strawberry.experimental.pydantic import type as pydantic_type
from strawberry.types.info import Info

from pyeed.ingest.model.annotation import Annotation
from pyeed.ingest.model.goannotation import GOAnnotation
from pyeed.ingest.model.molecule import Molecule
from pyeed.ingest.model.protein import Protein
from pyeed.ingest.model.reaction import Reaction
from pyeed.ingest.model.taxon import Taxon

from .context import GraphQLContext

if TYPE_CHECKING:
    from .types import ProteinType


@pydantic_type(Molecule)
class MoleculeType:
    id: strawberry.auto
    name: strawberry.auto
    smiles: strawberry.auto
    inchi: strawberry.auto


@pydantic_type(Reaction)
class ReactionType:
    id: strawberry.auto
    description: strawberry.auto
    reversible: strawberry.auto

    @strawberry.field(description="Substrates of the reaction")
    async def substrates(self, info: Info[GraphQLContext]) -> list[MoleculeType]:
        return await info.context.substrates_of_reaction.load(self.id)

    @strawberry.field(description="Products of the reaction")
    async def products(self, info: Info[GraphQLContext]) -> list[MoleculeType]:
        return await info.context.products_of_reaction.load(self.id)

    @strawberry.field(description="Proteins catalyzing the reaction")
    async def catalyzingProteins(
        self, info: Info[GraphQLContext]
    ) -> list[Annotated[ProteinType, strawberry.lazy(".types")]]:
        return await info.context.proteins_of_reaction.load(self.id)


@pydantic_type(Taxon)
class TaxonType:
    id: strawberry.auto
    scientific_name: strawberry.auto
    common_name: strawberry.auto
    rank: strawberry.auto
    hidden: strawberry.auto
    synonyms: strawberry.auto


@pydantic_type(Annotation)
class AnnotationType:
    id: strawberry.auto
    positions: strawberry.auto
    description: strawberry.auto


@pydantic_type(GOAnnotation)
class GOAnnotationType:
    id: strawberry.auto
    term: strawberry.auto
    definition: strawberry.auto


@pydantic_type(Protein)
class ProteinType:
    id: strawberry.auto
    sequence: strawberry.auto
    seq_length: strawberry.auto
    name: strawberry.auto
    mol_weight: strawberry.auto
    ec_numbers: strawberry.auto

    @strawberry.field(description="Vector of the protein")
    async def vector(self, info: Info[GraphQLContext]) -> list[float]:
        vector = await info.context.protein_vector_loader.load(key=self.id)
        return vector

    @strawberry.field(description="Organism the protein originates from")
    async def organism(self, info: Info[GraphQLContext]) -> TaxonType | None:
        return await info.context.organism_of_protein.load(self.id)

    @strawberry.field(description="Reactions catalyzed by the protein")
    async def reactions(self, info: Info[GraphQLContext]) -> list[ReactionType]:
        return await info.context.reactions_of_protein.load(self.id)

    @strawberry.field(description="GO annotations of the protein")
    async def goAnnotations(self, info: Info[GraphQLContext]) -> list[GOAnnotationType]:
        return await info.context.go_annotations_of_protein.load(self.id)

    @strawberry.field(description="Annotations of the protein (e.g. domains, families, etc.)")
    async def annotations(self, info: Info[GraphQLContext]) -> list[AnnotationType]:
        return await info.context.annotations_of_protein.load(self.id)
