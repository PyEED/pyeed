from dataclasses import dataclass

from neo4j import AsyncDriver

from pyeed.db.milvus import VectorDB


@dataclass
class GraphQLContext:
    """Injected context available in all resolvers."""

    neo4j_driver: AsyncDriver
    milvus_client: VectorDB
