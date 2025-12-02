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


@pydantic_type(Molecule)
class MoleculeType:
    inchi_key: strawberry.auto
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
        instances = await Reaction.get_related(
            Molecule,
            driver=info.context.neo4j_driver,
            id=self.id,
        )
        return [MoleculeType(**instance.model_dump()) for instance in instances]

    @strawberry.field(description="Products of the reaction")
    async def products(self, info: Info[GraphQLContext]) -> list[MoleculeType]:
        instances = await Reaction.get_related(
            Molecule,
            driver=info.context.neo4j_driver,
            id=self.id,
        )
        return [MoleculeType(**instance.model_dump()) for instance in instances]


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
    # annotation_type: strawberry.auto
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

    @strawberry.field(description="Organism the protein originates from")
    async def organism(self, info: Info[GraphQLContext]) -> TaxonType | None:
        instance = await Protein.get_related(
            Taxon,
            driver=info.context.neo4j_driver,
            id=self.id,
            direction="out",
        )
        return TaxonType(**instance[0].model_dump()) if instance else None

    @strawberry.field(description="Reactions catalyzed by the protein")
    async def reactions(self, info: Info[GraphQLContext]) -> list[ReactionType]:
        instances = await Protein.get_related(
            Reaction,
            driver=info.context.neo4j_driver,
            id=self.id,
        )
        return [ReactionType(**instance.model_dump()) for instance in instances]

    @strawberry.field(description="GO annotations of the protein")
    async def goAnnotations(self, info: Info[GraphQLContext]) -> list[GOAnnotationType]:
        instances = await Protein.get_related(
            GOAnnotation,
            driver=info.context.neo4j_driver,
            id=self.id,
        )
        return [GOAnnotationType(**instance.model_dump()) for instance in instances]

    @strawberry.field(description="Annotations of the protein (e.g. domains, families, etc.)")
    async def annotations(self, info: Info[GraphQLContext]) -> list[AnnotationType]:
        instances = await Protein.get_related(
            Annotation,
            driver=info.context.neo4j_driver,
            id=self.id,
        )
        return [AnnotationType(**instance.model_dump()) for instance in instances]
