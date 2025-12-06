from strawberry.asgi import GraphQL

from pyeed.db.milvus import get_async_milvus_client
from pyeed.db.neo4j import get_async_driver
from pyeed.graphql.loaders import (
    annotations_of_protein,
    get_protein_vector_loader,
    go_annotations_of_protein,
    organism_of_protein,
    products_of_reaction,
    proteins_of_reaction,
    reactions_of_protein,
    substrates_of_reaction,
)

from .context import GraphQLContext
from .schema import schema


class CustomGraphQL(GraphQL):
    """Custom GraphQL ASGI app with context injection."""

    async def get_context(self, request, response):
        """Override to provide custom context.

        Contains:
        - Neo4j driver
        - Milvus client
        - Protein vector loader
        - loaders to get multiple relationships at once for each db entity
        """
        neo4j_driver = get_async_driver()
        milvus_client = get_async_milvus_client()
        return GraphQLContext(
            neo4j_driver=neo4j_driver,
            milvus_client=milvus_client,
            # Protein relationships
            protein_vector_loader=get_protein_vector_loader(milvus_client=milvus_client),
            organism_of_protein=organism_of_protein(neo4j_driver=neo4j_driver),
            reactions_of_protein=reactions_of_protein(neo4j_driver=neo4j_driver),
            go_annotations_of_protein=go_annotations_of_protein(neo4j_driver=neo4j_driver),
            annotations_of_protein=annotations_of_protein(neo4j_driver=neo4j_driver),
            # Reaction relationships
            proteins_of_reaction=proteins_of_reaction(neo4j_driver=neo4j_driver),
            substrates_of_reaction=substrates_of_reaction(neo4j_driver=neo4j_driver),
            products_of_reaction=products_of_reaction(neo4j_driver=neo4j_driver),
        )


app = CustomGraphQL(schema)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("pyeed.graphql.app:app", host="0.0.0.0", port=8123, reload=True)
