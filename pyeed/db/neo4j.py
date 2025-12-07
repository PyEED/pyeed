from __future__ import annotations

import asyncio
import os
from collections import defaultdict
from dataclasses import dataclass

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


@dataclass
class Relationship:
    from_label: str
    rel_type: str
    to_label: str


@dataclass
class Property:
    name: str
    is_index: bool = False


@dataclass
class Label:
    name: str
    properties: list[Property]


@dataclass
class Schema:
    labels: list[Label]
    relationships: list[Relationship]


async def get_labels(driver: AsyncDriver) -> list[Label]:
    """
    Get the labels of the database.

    Args:
        driver: Neo4j async driver.

    Returns:
        list[Label]: The labels of the database.
    """
    query = """
    CALL apoc.meta.schema() YIELD value AS schemaMap
    UNWIND keys(schemaMap) AS label
    WITH label, schemaMap[label] AS data
    WHERE data.type = 'node'
    UNWIND keys(data.properties) AS prop
    WITH label, prop, data.properties[prop] AS propData
    RETURN
      label              AS label_name,
      prop               AS property_name,
      coalesce(propData.indexed, false) AS is_index
    """
    async with driver.session() as session:

        async def _tx(tx: AsyncManagedTransaction):
            result = await tx.run(query)
            return [record async for record in result]

        rows = await session.execute_read(_tx)

    by_label: defaultdict[str, list[Property]] = defaultdict(list)

    for row in rows:
        by_label[row["label_name"]].append(
            Property(
                name=row["property_name"],
                is_index=row["is_index"],
            )
        )

    return [Label(name=label_name, properties=props) for label_name, props in by_label.items()]


async def get_relationships(driver: AsyncDriver) -> list[Relationship]:
    """
    Get the relationships of the database.

    Args:
        driver: Neo4j async driver.

    Returns:
        list[Relationship]: The relationships of the database.
    """
    query = """
    CALL apoc.meta.graph()
    YIELD relationships
    UNWIND relationships AS rel
    UNWIND labels(startNode(rel)) AS from_label
    UNWIND labels(endNode(rel))   AS to_label
    WITH DISTINCT from_label, type(rel) AS rel_type, to_label
    RETURN from_label, rel_type, to_label
    ORDER BY from_label, rel_type, to_label;
    """

    async with driver.session() as session:

        async def _tx(tx: AsyncManagedTransaction):
            result = await tx.run(query)
            return [record async for record in result]

        rows = await session.execute_read(_tx)

    relationships = []
    for row in rows:
        relationships.append(
            Relationship(
                from_label=row.get("from_label"),
                rel_type=row.get("rel_type"),
                to_label=row.get("to_label"),
            )
        )

    return relationships


async def get_schema(driver: AsyncDriver) -> Schema:
    """
    Get the schema of the database.

    Args:
        driver: Neo4j async driver.

    Returns:
        Schema: The schema of the database.
    """
    labels = await get_labels(driver)
    relationships = await get_relationships(driver)
    return Schema(labels=labels, relationships=relationships)


if __name__ == "__main__":
    from rich import print

    driver = get_async_driver()

    async def main() -> None:
        schema = await get_schema(driver)
        print(schema)
        await driver.close()

    asyncio.run(main())
