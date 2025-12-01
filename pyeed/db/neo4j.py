from __future__ import annotations

import asyncio
import os

import dotenv
from loguru import logger
from neo4j import AsyncDriver, AsyncGraphDatabase, AsyncManagedTransaction, Driver, GraphDatabase

from ..ingest.model.pyeedbase import BaseNode
from ..ingest.model.utils import collect_schema


def get_async_driver(
    uri: str | None = None, user: str | None = None, password: str | None = None
) -> AsyncDriver:
    dotenv.load_dotenv()
    uri = uri or os.getenv("NEO4J_URI")
    user = user or os.getenv("NEO4J_USER")
    password = password or os.getenv("NEO4J_PASSWORD")
    if not (uri and user and password):
        raise ValueError(
            "URI, user, and password must be provided or set in env "
            f"(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD): {uri}, {user}, {password}"
        )
    return AsyncGraphDatabase.driver(uri, auth=(user, password))


def get_driver(
    uri: str | None = None, user: str | None = None, password: str | None = None
) -> Driver:
    dotenv.load_dotenv()
    uri = uri or os.getenv("NEO4J_URI")
    user = user or os.getenv("NEO4J_USER")
    password = password or os.getenv("NEO4J_PASSWORD")
    if not (uri and user and password):
        raise ValueError(
            "URI, user, and password must be provided or set in env "
            f"(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD): {uri}, {user}, {password}"
        )
    return GraphDatabase.driver(uri, auth=(user, password))


async def sync_schema(models: list[type[BaseNode]], driver: AsyncDriver) -> None:
    """
    Synchronize database schema from model definitions.

    Args:
        models: List of model classes to synchronize.
        driver: Neo4j async driver.

    Returns:
        None
    """
    logger.info("Synchronizing schema from model definitions...")
    added_indices = []
    async with driver.session() as session:
        for spec in collect_schema(models):
            cypher = (
                f"CREATE CONSTRAINT `uniq_{spec.label}_{spec.prop}` IF NOT EXISTS "
                f"FOR (n:`{spec.label}`) REQUIRE n.`{spec.prop}` IS UNIQUE"
            )
            added_indices.append(spec)

            async def _tx(tx: AsyncManagedTransaction, query: str = cypher) -> None:
                await tx.run(query)

            await session.execute_write(_tx)

    logger.info(f"Synchronized schema for the following indices: {added_indices}")

    # async def async_query_iter(self, query: str, **params: Any) -> AsyncIterator[dict[str, Any]]:
    #     """
    #     Iterate over the results of a query.

    #     Args:
    #         query: The query to execute.
    #         params: The parameters to pass to the query.

    #     Returns:
    #         An iterator over the results of the query.
    #     """
    #     async with self.async_driver.session() as session:
    #         result = await session.run(query, **params)
    #         async for record in result:
    #             yield record.data()

    # async def async_value_iter(self, query: str, key: str, **params: Any) -> AsyncIterator[Any]:
    #     """
    #     Iterate over the values of a key from the results of a query.

    #     Args:
    #         query: The query to execute.
    #         key: The key to iterate over.
    #         params: The parameters to pass to the query.

    #     Returns:
    #         An iterator over the values of the key.
    #     """
    #     async with self.async_driver.session() as session:
    #         result = await session.run(query, **params)
    #         async for record in result:
    #             yield record.value(key)


async def main() -> None:
    uri = "neo4j://127.0.0.1:7687"
    user = "neo4j"
    password = "12345678"
    driver = get_async_driver(uri, user, password)
    try:
        await sync_schema(MODEL_CLASSES, driver)
    finally:
        await driver.close()
    print("Done")


if __name__ == "__main__":
    from pyeed.ingest.model import MODEL_CLASSES

    asyncio.run(main())
