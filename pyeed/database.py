from __future__ import annotations

import logging
import os
import re
from collections.abc import AsyncIterator, Iterable
from typing import Any

import dotenv
from neo4j import AsyncGraphDatabase, AsyncSession, GraphDatabase

from .model.pyeedbase import LabelProperty, PyeedBase
from .model.utility import collect_schema

logger = logging.getLogger(__name__)


class Database:
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
                "(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)"
            )
        self.async_driver = AsyncGraphDatabase.driver(uri, auth=(user, password))
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
        self._vec_index_cache: set[str] = set()

    async def close(self) -> None:
        await self.async_driver.close()

    def verify_connection(self) -> None:
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
        with self.driver.session() as session:
            result = session.run(query, **params)
            return [record.data() for record in result]

    async def async_query_iter(self, query: str, **params: Any) -> AsyncIterator[dict[str, Any]]:
        """
        Iterate over the results of a query.
        """
        async with self.async_driver.session() as session:
            result = await session.run(query, **params)
            async for record in result:
                yield record.data()

    async def async_value_iter(self, query: str, key: str, **params: Any) -> AsyncIterator[Any]:
        """
        Iterate over the values of a key from the results of a query.
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

    async def _ensure_vector_index(
        self,
        session: AsyncSession,
        prop: str,
        dims: int,
        sim: str = "cosine",
    ) -> None:
        """
        Ensure a vector index exists on Embedding.`prop`.

        Notes:
            - Neo4j does NOT allow parameterizing index names or property identifiers.
            - We inline both and only parameterize values (dims, sim).
        """
        logger.info(f"Setting up vector index for {prop} with dimensions {dims}...")
        if prop in self._vec_index_cache:
            return

        # sanitize for safety
        safe_prop = prop.replace("`", "``")  # escape backticks in property
        # keep index name simple (letters, digits, underscores)
        idx_name = f"vx_Embedding__{prop}"
        idx_name = re.sub(r"[^A-Za-z0-9_]", "_", idx_name)

        q = (
            f"CREATE VECTOR INDEX {idx_name} IF NOT EXISTS "
            f"FOR (n:Embedding) ON (n.`{safe_prop}`) "
            "OPTIONS {indexConfig: { `vector.dimensions`: $dim, "
            "`vector.similarity_function`: $sim }}"
        )
        await session.run(q, dim=dims, sim=sim)
        self._vec_index_cache.add(prop)

    async def bulk_upsert(
        self,
        nodes: list[dict[str, Any]],
        edges: list[dict[str, Any]],
        tx_size: int = 5000,
    ) -> None:
        """
        Insert or update a batch of nodes and relationships.

        Args:
            nodes: List of node dicts, each with keys:
                - "label": Node label (str).
                - "key": (unique_field, unique_value).
                - "props": Dict of node properties.
            edges: List of edge dicts, each with keys:
                - "type": Relationship type (str).
                - "src": (label, key_field, key_value) of source node.
                - "dst": (label, key_field, key_value) of target node.
            tx_size: Maximum number of rows per transaction batch.

        Notes:
            - Nodes are grouped by label and upserted with MERGE.
            - Relationships are grouped by type and endpoints and merged.
            - Vector indexes for embeddings are created automatically if needed.
        """
        logger.info(f"Bulk upserting {len(nodes)} nodes and {len(edges)} edges...")
        async with self.async_driver.session() as session:
            # 0) Prepare vector indexes (Embedding vec__* props)
            vec_props: dict[str, int] = {}
            for n in nodes:
                if n.get("label") != "Embedding":
                    continue
                props = n.get("props", {})
                dims = props.get("n_dims")
                if isinstance(dims, int):
                    for k in props:
                        if k.startswith("vec__"):
                            vec_props[k] = dims
            for prop, dim in vec_props.items():
                await self._ensure_vector_index(session, prop, dim)

            # Nodes by label
            buckets: dict[str, tuple[str, list[dict[str, Any]]]] = {}
            for n in nodes:
                lbl = n["label"]
                key_name, key_val = n["key"]
                buckets.setdefault(lbl, (key_name, []))[1].append(
                    {"key": key_val, "props": n["props"]}
                )

            for lbl, (kname, rows) in buckets.items():
                for i in range(0, len(rows), tx_size):
                    chunk = rows[i : i + tx_size]
                    q = f"""
                    UNWIND $rows AS r
                    MERGE (n:`{lbl}` {{ `{kname}`: r.key }})
                    SET n += r.props
                    """
                    await session.run(q, rows=chunk)

            # 2) Relationships grouped by (etype, sl, sk, dl, dk)
            groups: dict[tuple[str, str, str, str, str], list[dict[str, Any]]] = {}
            for e in edges:
                sl, sk, sv = e["src"]
                dl, dk, dv = e["dst"]
                key = (e["type"], sl, sk, dl, dk)
                groups.setdefault(key, []).append({"sv": sv, "dv": dv})

            for (etype, sl, sk, dl, dk), rows in groups.items():
                for i in range(0, len(rows), tx_size):
                    chunk = rows[i : i + tx_size]
                    q = f"""
                    UNWIND $rows AS r
                    MATCH (s:`{sl}` {{ `{sk}`: r.sv }})
                    MATCH (d:`{dl}` {{ `{dk}`: r.dv }})
                    MERGE (s)-[:`{etype}`]->(d)
                    """
                    await session.run(q, rows=chunk)

    async def save(self, node: PyeedBase) -> None:
        """
        Insert or update a single root node and its entire subtree.

        Args:
            node: A BaseNode instance (e.g. Protein) with nested child nodes.
        """
        logger.info(f"Upserting node {node.__class__.__name__}...")
        nodes, edges = node.graphify()
        await self.bulk_upsert(nodes, edges)

    async def save_many(
        self,
        roots: Iterable[PyeedBase],
        tx_size: int = 5000,
    ) -> None:
        """
        Insert or update multiple root nodes and their subtrees in one call.

        Args:
            roots: Iterable of PyeedBase instances to upsert.
            tx_size: Maximum number of rows per transaction batch.
        """
        logger.info("Upserting nodes...")
        all_nodes: list[dict[str, Any]] = []
        all_edges: list[dict[str, Any]] = []
        for r in roots:
            n, e = r.graphify()
            all_nodes.extend(n)
            all_edges.extend(e)
        await self.bulk_upsert(all_nodes, all_edges, tx_size=tx_size)

    async def attach(
        self,
        parent: type[PyeedBase],
        parents_to_children: dict[str | int, list[PyeedBase]],
        tx_size: int = 5000,
    ) -> None:
        """
        Attach child nodes to existing parents using the child's EDGES mapping.

        Args:
            parent: Parent node class (e.g., Protein).
            parents_to_children: Mapping of parent unique values to lists of children.
                Example: { "P01234": [Embedding(...), Annotation(...)] }
            tx_size: Maximum number of rows per transaction batch.

        Raises:
            ValueError: If a child node has no applicable EDGES rule for this parent,
                        or multiple ambiguous rules exist.

        Notes:
            - Each child is inserted along with its own subtree.
            - The relationship type is resolved from child_cls.EDGES via
              child_cls.resolve_edge(parent_label=<parent>, field_name=None).
              This requires a single unambiguous rule for this parent.
            - The parent unique field is resolved automatically from LabelProperty(unique=True).
        """
        logger.info(
            "Attaching children to %s: %d parents",
            parent.__name__,
            len(parents_to_children),
        )
        plabel = parent.__name__
        pkey = self._unique_key_of(parent)

        all_nodes: list[dict[str, Any]] = []
        all_edges: list[dict[str, Any]] = []

        for pval, children in parents_to_children.items():
            for child in children:
                child_cls = type(child)
                # resolve rel type from child's EDGES; field_name=None for attach()
                rel_type = child_cls.resolve_edge(parent_label=plabel, field_name=None)

                # child's own subtree
                cnodes, cedges = child.graphify()
                all_nodes.extend(cnodes)
                all_edges.extend(cedges)

                # explicit parent -> child edge
                ckey = child.get_unique_model_field()
                cval = getattr(child, ckey)
                all_edges.append(
                    {
                        "type": rel_type,
                        "src": (plabel, pkey, pval),
                        "dst": (child_cls.__name__, ckey, cval),
                    }
                )

        await self.bulk_upsert(all_nodes, all_edges, tx_size=tx_size)

    # Helpers
    @staticmethod
    def _unique_key_of(model_cls: type[PyeedBase]) -> str:
        """
        Get the unique key field of a model class.

        Args:
            model_cls: A PyeedBase subclass.

        Returns:
            The name of the field marked with `NodeHint(unique=True)`.

        Raises:
            ValueError: If the model class has no unique field.
        """
        logger.debug(f"Finding unique key for {model_cls.__name__}...")
        for fname, finfo in model_cls.model_fields.items():
            meta = getattr(finfo, "metadata", []) or ()
            for m in meta:
                if isinstance(m, LabelProperty) and m.unique:
                    return str(fname)
        raise ValueError(f"No unique field declared for {model_cls.__name__}")
