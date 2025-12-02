from strawberry.asgi import GraphQL

from pyeed.db.milvus import get_async_milvus_client
from pyeed.db.neo4j import get_async_driver

from .context import GraphQLContext
from .schema import schema


class CustomGraphQL(GraphQL):
    """Custom GraphQL ASGI app with context injection."""

    async def get_context(self, request, response):
        """Override to provide custom context."""
        return GraphQLContext(
            neo4j_driver=get_async_driver(),
            milvus_client=get_async_milvus_client(),
        )


app = CustomGraphQL(schema)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("pyeed.graphql.app:app", host="0.0.0.0", port=8123, reload=True)
