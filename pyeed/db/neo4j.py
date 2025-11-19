from __future__ import annotations

import os
from collections.abc import AsyncIterator
from typing import Any

import dotenv
from loguru import logger
from neo4j import AsyncGraphDatabase, GraphDatabase
from pandas.core.common import defaultdict

from ..ingest.model.pyeedbase import PyeedBase
from ..ingest.model.utils import collect_schema


class GraphDB:
    def __init__(
        self,
        uri: str | None = None,
        user: str | None = None,
        password: str | None = None,
    ) -> None:
        dotenv.load_dotenv()
        uri = uri or os.getenv("NEO4J_URI")
        user = user or os.getenv("NEO4J_USER")
        password = password or os.getenv("NEO4J_PASSWORD")
        if not (uri and user and password):
            raise ValueError(
                "URI, user, and password must be provided or set in env "
                f"(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD): {uri}, {user}, {password}"
            )
        self.uri = uri
        self.async_driver = AsyncGraphDatabase.driver(uri, auth=(user, password))
        self.driver = GraphDatabase.driver(uri, auth=(user, password))

    async def close(self) -> None:
        """
        Close the database connection.
        """
        await self.async_driver.close()

    def verify_connection(self) -> None:
        """
        Verify the database connection.
        Prints "ok" if the connection is successful.
        """
        with self.driver.session() as session:
            result = session.run("RETURN 1 AS ok")
            record = result.single()
            if not record or record["ok"] != 1:
                raise ConnectionError("Failed to verify Neo4j connection")
            print("ok")

    def query(
        self,
        query: str,
        **params: Any,
    ) -> list[dict[str, Any]]:
        """
        Query the database.

        Args:
            query: The query to execute.
            params: The parameters to pass to the query.

        Returns:
            The results of the query.
        """
        with self.driver.session() as session:
            return session.run(query, **params).data()

    async def async_query_iter(self, query: str, **params: Any) -> AsyncIterator[dict[str, Any]]:
        """
        Iterate over the results of a query.

        Args:
            query: The query to execute.
            params: The parameters to pass to the query.

        Returns:
            An iterator over the results of the query.
        """
        async with self.async_driver.session() as session:
            result = await session.run(query, **params)
            async for record in result:
                yield record.data()

    async def async_value_iter(self, query: str, key: str, **params: Any) -> AsyncIterator[Any]:
        """
        Iterate over the values of a key from the results of a query.

        Args:
            query: The query to execute.
            key: The key to iterate over.
            params: The parameters to pass to the query.

        Returns:
            An iterator over the values of the key.
        """
        async with self.async_driver.session() as session:
            result = await session.run(query, **params)
            async for record in result:
                yield record.value(key)

    async def sync_schema(self, models: list[type[PyeedBase]]) -> None:
        """
        Synchronize database schema from model definitions.

        Creates unique constraints and btree indexes according to `NodeHint` metadata.
        This is idempotent and can be run multiple times without side effects.

        Args:
            models: Iterable of BaseNode subclasses to inspect for schema metadata.
        """
        logger.info("Creating unique constraints ...")
        uniques, btrees, _ = collect_schema(models)
        async with self.async_driver.session() as session:
            # Unique constraints
            for uniq in uniques:
                q = (
                    f"CREATE CONSTRAINT `uniq_{uniq.label}_{uniq.prop}` IF NOT EXISTS "
                    f"FOR (n:`{uniq.label}`) REQUIRE n.`{uniq.prop}` IS UNIQUE"
                )
                logger.debug("Schema: %s", q)
                await session.run(q)

            # B-tree indexes
            for btree in btrees:
                q = (
                    f"CREATE INDEX `idx_{btree.label}_{btree.prop}` IF NOT EXISTS "
                    f"FOR (n:`{btree.label}`) ON (n.`{btree.prop}`)"
                )
                logger.debug("Schema: %s", q)
                await session.run(q)

    async def upsert_nodes(
        self,
        nodes: list[PyeedBase],
        tx_size: int = 5000,
    ) -> defaultdict[str, set[str]]:
        """Upsert PyeedBase nodes without edges.

        Groups nodes by label, uses MERGE on unique field with SET for all properties.

        Args:
            nodes: Pre-batched list of PyeedBase instances
            tx_size: Max rows per Neo4j UNWIND transaction

        Returns:
            List of (label, unique_id) for each processed node
        """

        tracking: defaultdict[str, set[str]] = defaultdict(set)

        if not nodes:
            return tracking

        # Group by label (class name)
        by_label: dict[str, list[PyeedBase]] = {}
        for node in nodes:
            label = type(node).__name__
            by_label.setdefault(label, []).append(node)

        async with self.async_driver.session() as session:
            for label, node_list in by_label.items():
                # Get unique field from first node
                unique_field = node_list[0].get_unique_model_field()

                # Build rows for UNWIND
                rows = []
                for node in node_list:
                    unique_value = str(getattr(node, unique_field))
                    props = node.to_dict()
                    rows.append({"key": unique_value, "props": props})

                # Batch by tx_size
                for i in range(0, len(rows), tx_size):
                    chunk = rows[i : i + tx_size]
                    query = f"""
                    UNWIND $rows AS r
                    MERGE (n:`{label}` {{ `{unique_field}`: r.key }})
                    SET n += r.props
                    """
                    await session.run(query, rows=chunk)

                    # simulate longer processing time to test queue backpressure
                    tracking[label].update(set(r["key"] for r in chunk))

        return tracking
