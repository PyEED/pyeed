from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Literal, Self, overload

from neo4j import AsyncDriver
from pydantic import BaseModel, ConfigDict, model_validator

from ...db.query_utils import (
    execute_read,
    execute_write,
    process_multiple_records,
    process_single_record,
)


@dataclass(frozen=True)
class LabelProperty:
    """Influences the Neo4j schema creation.

    Attributes:
        index: Whether this field should be indexed in Neo4j.
    """

    index: bool = False


class BaseNode(BaseModel):
    """Base class for all nodes in the Database.

    To mark a field as indexed, use ``Annotated`` with ``LabelProperty(index=True)``
    as a type hint annotation. For example::

        from typing import Annotated
        from pydantic import Field

        class MyNode(BaseNode):
            id: Annotated[str, LabelProperty(index=True)] = Field(
                ..., description="ID field"
            )

    Note:
        Only one field per class can be marked with ``index=True``.
    """

    model_config = ConfigDict(
        frozen=False,
        validate_assignment=True,
        use_enum_values=True,
        populate_by_name=True,
    )

    @model_validator(mode="after")
    def _check_indexed_fields(self) -> Self:
        """Validate that only one field is marked as indexed.

        Returns:
            Self: The validated instance.

        Raises:
            ValueError: If more than one field is marked with index=True.
        """
        cls = type(self)

        indexed = []
        for name, field in cls.model_fields.items():
            for m in field.metadata:
                if isinstance(m, LabelProperty) and m.index:
                    indexed.append(name)
                    break

        if len(indexed) > 1:
            raise ValueError(
                f"Only one field can be marked with LabelProperty(index=True). "
                f"Found {len(indexed)} indexed fields in {cls.__name__}: {indexed}"
            )

        return self

    # --------- Upsert methods --------- #

    async def _upsert[T: BaseNode](
        self: T,
        driver: AsyncDriver,
        session_kwargs: dict[str, Any] | None = None,
    ) -> None:
        """Upsert this node instance.

        Args:
            driver: Neo4j async driver.
            session_kwargs: Optional session keyword arguments.
        """
        await type(self)._bulk_upsert(
            driver,
            [self],
            session_kwargs=session_kwargs,
        )

    @classmethod
    async def _bulk_upsert[T: BaseNode](
        cls: type[T],
        driver: AsyncDriver,
        nodes: Iterable[T],
        session_kwargs: dict[str, Any] | None = None,
    ) -> None:
        """Upsert multiple nodes of this type in a single UNWIND batch.

        Args:
            driver: Neo4j async driver.
            nodes: Iterable of node instances to upsert.
            session_kwargs: Optional session keyword arguments.
        """
        rows = [node.model_dump() for node in nodes]
        if not rows:
            return

        index_field_name = cls.get_index_field()

        cypher = f"""
        UNWIND $rows AS row
        MERGE (n:`{cls.__name__}` {{ {index_field_name}: row.{index_field_name} }})
        SET n += row
        """

        await execute_write(
            driver=driver,
            query=cypher,
            rows=rows,
            session_kwargs=session_kwargs,
        )

    # --------- Relationship methods --------- #

    async def _relate[S: BaseNode, E: BaseNode](
        self: S,
        driver: AsyncDriver,
        rel_type: str,
        targets: Iterable[E] | E,
        *,
        session_kwargs: dict[str, Any] | None = None,
        rel_props: dict[str, Any] | None = None,
    ) -> None:
        """Create (self)-[rel_type]->(targets) relationships.

        Wrapper around `_bulk_create_relationships` for one start node.

        Args:
            driver: Neo4j async driver.
            rel_type: Relationship type.
            targets: Single target node or iterable of target nodes.
            session_kwargs: Optional session keyword arguments.
            rel_props: Optional relationship properties.
        """
        pairs = [(self, targets)] if isinstance(targets, BaseNode) else [(self, t) for t in targets]

        await type(self)._bulk_relate(
            driver=driver,
            rel_type=rel_type,
            pairs=pairs,
            session_kwargs=session_kwargs,
            rel_props=rel_props,
        )

    @classmethod
    async def _bulk_relate[S: BaseNode, E: BaseNode](
        cls: type[S],
        driver: AsyncDriver,
        rel_type: str,
        pairs: Iterable[tuple[S, E]],
        *,
        session_kwargs: dict[str, Any] | None = None,
        rel_props: dict[str, Any] | None = None,
    ) -> None:
        """Create (start:cls)-[rel_type]->(target) relationships for many pairs.

        Args:
            driver: Neo4j async driver.
            rel_type: Relationship type.
            pairs: Iterable of (start, target) pairs.
            session_kwargs: Optional session keyword arguments.
            rel_props: Optional relationship properties.

        Raises:
            ValueError: If pairs is empty, contains mixed types, or start nodes
                don't match the calling class.
        """
        pairs_list = list(pairs)
        if not pairs_list:
            return

        first_start, first_target = pairs_list[0]

        # Validate that start nodes match the calling class
        if not isinstance(first_start, cls):
            raise ValueError(
                f"_bulk_create_relationships must be called on the start class; "
                f"expected {cls.__name__}, got {type(first_start).__name__}"
            )

        # Get target class from first pair
        target_cls = type(first_target)

        # Validate all pairs have consistent types
        for i, (start, target) in enumerate(pairs_list, start=1):
            if not isinstance(start, cls):
                raise ValueError(
                    f"Inconsistent start node type at index {i}: "
                    f"expected {cls.__name__}, got {type(start).__name__}"
                )
            if not isinstance(target, target_cls):
                raise ValueError(
                    f"Inconsistent target node type at index {i}: "
                    f"expected {target_cls.__name__}, got {type(target).__name__}"
                )

        start_key = cls.get_index_field()
        target_label = target_cls.__name__
        target_key = target_cls.get_index_field()

        props = rel_props or {}
        rows = [
            {
                "from_id": getattr(start, start_key),
                "to_id": getattr(target, target_key),
                "props": props,
            }
            for start, target in pairs_list
        ]

        cypher = f"""
        UNWIND $rows AS row
        MATCH (a:`{cls.__name__}` {{ `{start_key}`: row.from_id }})
        MATCH (b:`{target_label}` {{ `{target_key}`: row.to_id }})
        MERGE (a)-[r:`{rel_type}`]->(b)
        SET r += row.props
        """

        await execute_write(
            driver=driver,
            query=cypher,
            rows=rows,
            session_kwargs=session_kwargs,
        )

    @classmethod
    def get_index_field(cls) -> str:
        """Returns the name of the field marked with LabelProperty(index=True).

        Returns:
            The name of the indexed field.

        Raises:
            ValueError: If no field is marked as indexed.
        """
        for field_name, field_info in cls.model_fields.items():
            for m in field_info.metadata:
                if isinstance(m, LabelProperty) and m.index:
                    return field_name

        raise ValueError(
            f"No field with LabelProperty(index=True) found in {cls.__name__}. "
            f"To mark a field as indexed, use: "
            f"Annotated[<type>, LabelProperty(index=True)] as the type hint "
            f"annotation."
        )

    # --------- Query methods --------- #

    @classmethod
    @overload
    async def get(
        cls,
        driver: AsyncDriver,
        *,
        id: str,
    ) -> Self | None:
        """Get a single node by its indexed field value."""
        ...

    @classmethod
    @overload
    async def get(
        cls,
        driver: AsyncDriver,
        *,
        ids: list[str],
    ) -> list[Self]:
        """Get multiple nodes by their indexed field values."""
        ...

    @classmethod
    async def get(
        cls,
        driver: AsyncDriver,
        *,
        id: str | None = None,
        ids: list[str] | None = None,
    ) -> Self | list[Self] | None:
        """Get node(s) by indexed field.

        Args:
            driver: Neo4j async driver.
            id: Single ID to fetch (returns one node or None).
            ids: Multiple IDs to fetch (returns list).

        Returns:
            Single node, list of nodes, or None.

        Raises:
            ValueError: If neither id nor ids is provided, or both are provided.
        """
        if id is not None and ids is not None:
            raise ValueError("Provide either 'id' or 'ids', not both")
        if id is None and ids is None:
            raise ValueError("Provide either 'id' or 'ids'")

        index_field = cls.get_index_field()
        label = cls.__name__

        if id is not None:
            # Single lookup
            query = f"MATCH (n:{label} {{{index_field}: $value}}) RETURN properties(n) AS node"
            data = await execute_read(
                driver=driver,
                query=query,
                params={"value": id},
                processor=process_single_record,
            )
            if not data:
                return None
            return cls(**data["node"])

        else:
            # Batch lookup
            query = (
                f"MATCH (n:{label}) WHERE n.{index_field} IN $values RETURN properties(n) AS node"
            )
            records = await execute_read(
                driver=driver,
                query=query,
                params={"values": ids},
                processor=process_multiple_records,
            )
            return [cls(**r["node"]) for r in records]

    @classmethod
    async def get_by(
        cls,
        driver: AsyncDriver,
        *,
        ranges: dict[str, tuple[Any | None, Any | None]] | None = None,
        **filters: Any,
    ) -> list[Self]:
        """Get nodes by arbitrary field filters with optional range support.

        Validates that filter keys are valid model fields.

        Args:
            driver: Neo4j async driver.
            ranges: Optional dict mapping field names to (min, max) tuples.
                None values in tuples are ignored. Example:
                {"seq_length": (100, 200), "mol_weight": (50000.0, None)}
            **filters: Field name to value mappings for exact matches.

        Returns:
            List of matching nodes.

        Raises:
            ValueError: If a filter key is not a valid field name, or if both
                filters and ranges are empty.

        Example:
            # Exact match
            proteins = await Protein.get_by(driver, name="Hemoglobin")

            # Range query
            proteins = await Protein.get_by(
                driver,
                ranges={"seq_length": (100, 200), "mol_weight": (50000.0, None)}
            )

            # Combined exact match and range
            proteins = await Protein.get_by(
                driver,
                name="Hemoglobin",
                ranges={"seq_length": (100, None)}
            )
        """
        valid_fields = set(cls.model_fields.keys())

        # Validate exact match filter keys
        invalid_keys = set(filters.keys()) - valid_fields
        if invalid_keys:
            raise ValueError(
                f"Invalid filter keys for {cls.__name__}: {invalid_keys}. "
                f"Valid fields: {valid_fields}"
            )

        # Validate range filter keys
        if ranges:
            invalid_range_keys = set(ranges.keys()) - valid_fields
            if invalid_range_keys:
                raise ValueError(
                    f"Invalid range keys for {cls.__name__}: {invalid_range_keys}. "
                    f"Valid fields: {valid_fields}"
                )

        if not filters and not ranges:
            raise ValueError("At least one filter or range must be provided")

        label = cls.__name__
        conditions: list[str] = []
        params: dict[str, Any] = {}

        # Build exact match conditions
        for key, value in filters.items():
            conditions.append(f"n.`{key}` = ${key}")
            params[key] = value

        # Build range conditions
        if ranges:
            for field, (min_val, max_val) in ranges.items():
                if min_val is not None:
                    param_name = f"{field}_min"
                    conditions.append(f"n.`{field}` >= ${param_name}")
                    params[param_name] = min_val
                if max_val is not None:
                    param_name = f"{field}_max"
                    conditions.append(f"n.`{field}` <= ${param_name}")
                    params[param_name] = max_val

        where_clause = " AND ".join(conditions)
        query = f"MATCH (n:`{label}`) WHERE {where_clause} RETURN properties(n) AS node"

        records = await execute_read(
            driver=driver,
            query=query,
            params=params,
            processor=process_multiple_records,
        )
        return [cls(**r["node"]) for r in records]

    @classmethod
    async def get_all(
        cls,
        driver: AsyncDriver,
        *,
        limit: int = 100,
        offset: int = 0,
    ) -> list[Self]:
        """Get all nodes of this type with pagination.

        Args:
            driver: Neo4j async driver.
            limit: Maximum number of results.
            offset: Number of results to skip.

        Returns:
            List of nodes.
        """
        label = cls.__name__
        query = f"MATCH (n:{label}) RETURN properties(n) AS node SKIP $offset LIMIT $limit"

        records = await execute_read(
            driver=driver,
            query=query,
            params={"limit": limit, "offset": offset},
            processor=process_multiple_records,
        )
        return [cls(**r["node"]) for r in records]

    @classmethod
    async def count(cls, driver: AsyncDriver) -> int:
        """Count all nodes of this type.

        Args:
            driver: Neo4j async driver.

        Returns:
            Number of nodes.
        """
        label = cls.__name__
        query = f"MATCH (n:{label}) RETURN count(n) AS count"

        data = await execute_read(
            driver=driver,
            query=query,
            params={},
            processor=process_single_record,
        )
        return data["count"] if data else 0

    @overload
    @classmethod
    async def get_related[S: BaseNode, E: BaseNode](
        cls: type[S],
        target_cls: type[E],
        *,
        driver: AsyncDriver,
        id: str,
        rel_type: str | None = None,
        direction: Literal["out", "in", "both"] = "out",
    ) -> list[E]:
        """Get related nodes for a single start node."""
        ...

    @overload
    @classmethod
    async def get_related[S: BaseNode, E: BaseNode](
        cls: type[S],
        target_cls: type[E],
        *,
        driver: AsyncDriver,
        ids: list[str],
        rel_type: str | None = None,
        direction: Literal["out", "in", "both"] = "out",
    ) -> dict[str, list[E]]:
        """Get related nodes for multiple start nodes.

        Returns:
            Dict mapping start node ID to list of related nodes.
        """
        ...

    @classmethod
    async def get_related[S: BaseNode, E: BaseNode](
        cls: type[S],
        target_cls: type[E],
        *,
        driver: AsyncDriver,
        id: str | None = None,
        ids: list[str] | None = None,
        rel_type: str | None = None,
        direction: Literal["out", "in", "both"] = "out",
    ) -> list[E] | dict[str, list[E]]:
        """Get one-hop related nodes of another class.

        Args:
            target_cls: Target node class (must inherit from BaseNode).
            driver: Neo4j async driver.
            id: Single start node ID (returns list of related nodes).
            ids: Multiple start node IDs (returns dict mapping ID to list of related nodes).
            rel_type: Optional relationship type to match. If None, any type is matched.
            direction: Relationship direction: 'out', 'in', or 'both'.

        Returns:
            If `id` provided: List of target_cls instances.
            If `ids` provided: Dict mapping start node ID to list of target_cls instances.

        Raises:
            ValueError: If neither id nor ids is provided, or both are provided.
        """
        if id is not None and ids is not None:
            raise ValueError("Provide either 'id' or 'ids', not both")
        if id is None and ids is None:
            raise ValueError("Provide either 'id' or 'ids'")

        start_label = cls.__name__
        start_key = cls.get_index_field()
        target_label = target_cls.__name__

        dir_normalized = direction.lower()
        if dir_normalized == "out":
            rel_pattern = f"-[r:`{rel_type}`]->" if rel_type else "-[r]->"
        elif dir_normalized == "in":
            rel_pattern = f"<-[r:`{rel_type}`]-" if rel_type else "<-[r]-"
        elif dir_normalized == "both":
            rel_pattern = f"-[r:`{rel_type}`]-" if rel_type else "-[r]-"
        else:
            raise ValueError("direction must be 'out', 'in', or 'both'")

        if id is not None:
            # Single lookup
            query = f"""
            MATCH (start:`{start_label}` {{ `{start_key}`: $id }}){rel_pattern}(target:`{target_label}`)
            RETURN properties(target) AS node
            """
            records = await execute_read(
                driver=driver,
                query=query,
                params={"id": id},
                processor=process_multiple_records,
            )
            return [target_cls(**r["node"]) for r in records]

        else:
            # Batch lookup
            query = f"""
            MATCH (start:`{start_label}`){rel_pattern}(target:`{target_label}`)
            WHERE start.`{start_key}` IN $ids
            RETURN start.`{start_key}` AS start_id, properties(target) AS node
            """
            records = await execute_read(
                driver=driver,
                query=query,
                params={"ids": ids},
                processor=process_multiple_records,
            )
            # Group by start_id
            result: dict[str, list[E]] = {}
            for r in records:
                start_id = r["start_id"]
                if start_id not in result:
                    result[start_id] = []
                result[start_id].append(target_cls(**r["node"]))
            return result

    @classmethod
    @overload
    async def get_filtered(
        cls,
        driver: AsyncDriver,
        *,
        id: str,
        ranges: dict[str, tuple[Any | None, Any | None]] | None = None,
        **filters: Any,
    ) -> Self | None:
        """Get a single node by index field with additional filters and ranges."""
        ...

    @classmethod
    @overload
    async def get_filtered(
        cls,
        driver: AsyncDriver,
        *,
        ids: list[str],
        ranges: dict[str, tuple[Any | None, Any | None]] | None = None,
        **filters: Any,
    ) -> list[Self]:
        """Get multiple nodes by index field with additional filters and ranges."""
        ...

    @classmethod
    @overload
    async def get_filtered(
        cls,
        driver: AsyncDriver,
        *,
        ranges: dict[str, tuple[Any | None, Any | None]] | None = None,
        **filters: Any,
    ) -> list[Self]:
        """Get nodes by filters and ranges only (no index field filtering)."""
        ...

    @classmethod
    async def get_filtered(
        cls,
        driver: AsyncDriver,
        *,
        id: str | None = None,
        ids: list[str] | None = None,
        ranges: dict[str, tuple[Any | None, Any | None]] | None = None,
        **filters: Any,
    ) -> Self | list[Self] | None:
        """Get node(s) by index field with additional filters and ranges.

        Combines index field filtering (from `get`) with arbitrary field filtering
        and ranges (from `get_by`). Can also be used without index field filtering
        (like `get_by`).

        Args:
            driver: Neo4j async driver.
            id: Single index field value to match (returns one node or None).
            ids: Multiple index field values to match (returns list).
            ranges: Optional dict mapping field names to (min, max) tuples.
                None values in tuples are ignored.
            **filters: Field name to value mappings for exact matches.

        Returns:
            Single node (if `id` provided), list of nodes (if `ids` provided or
            neither `id` nor `ids` provided), or None.

        Raises:
            ValueError: If both `id` and `ids` are provided, or if no filtering
                criteria is provided (at least one of: `id`, `ids`, `filters`, or
                `ranges` must be provided), or if a filter/range key is not a
                valid field name.

        Example:
            # Single ID with additional filter
            protein = await Protein.get_filtered(
                driver,
                id="P12345",
                name="Hemoglobin"
            )

            # Multiple IDs with range filter
            proteins = await Protein.get_filtered(
                driver,
                ids=["P12345", "P0CW62"],
                ranges={"seq_length": (100, 200)}
            )

            # Filter only (no index field) - like get_by
            proteins = await Protein.get_filtered(
                driver,
                name="Hemoglobin",
                ranges={"seq_length": (100, 200)}
            )

            # Range only (no index field)
            proteins = await Protein.get_filtered(
                driver,
                ranges={"mol_weight": (50000.0, None)}
            )
        """
        if id is not None and ids is not None:
            raise ValueError("Provide either 'id' or 'ids', not both")

        # Require at least one filtering criterion
        if id is None and ids is None and not filters and not ranges:
            raise ValueError(
                "At least one filtering criterion must be provided: id, ids, filters, or ranges"
            )

        valid_fields = set(cls.model_fields.keys())
        index_field = cls.get_index_field()
        label = cls.__name__

        # Validate exact match filter keys
        invalid_keys = set(filters.keys()) - valid_fields
        if invalid_keys:
            raise ValueError(
                f"Invalid filter keys for {cls.__name__}: {invalid_keys}. "
                f"Valid fields: {valid_fields}"
            )

        # Validate range filter keys
        if ranges:
            invalid_range_keys = set(ranges.keys()) - valid_fields
            if invalid_range_keys:
                raise ValueError(
                    f"Invalid range keys for {cls.__name__}: {invalid_range_keys}. "
                    f"Valid fields: {valid_fields}"
                )

        conditions: list[str] = []
        params: dict[str, Any] = {}

        # Build index field condition (only if id/ids provided)
        if id is not None:
            conditions.append(f"n.`{index_field}` = $index_value")
            params["index_value"] = id
        elif ids is not None:
            conditions.append(f"n.`{index_field}` IN $index_values")
            params["index_values"] = ids

        # Build exact match conditions
        for key, value in filters.items():
            conditions.append(f"n.`{key}` = ${key}")
            params[key] = value

        # Build range conditions
        if ranges:
            for field, (min_val, max_val) in ranges.items():
                if min_val is not None:
                    param_name = f"{field}_min"
                    conditions.append(f"n.`{field}` >= ${param_name}")
                    params[param_name] = min_val
                if max_val is not None:
                    param_name = f"{field}_max"
                    conditions.append(f"n.`{field}` <= ${param_name}")
                    params[param_name] = max_val

        where_clause = " AND ".join(conditions)
        query = f"MATCH (n:`{label}`) WHERE {where_clause} RETURN properties(n) AS node"

        if id is not None:
            # Single lookup - return single node or None
            records = await execute_read(
                driver=driver,
                query=query,
                params=params,
                processor=process_multiple_records,
            )
            if not records:
                return None
            return cls(**records[0]["node"])
        else:
            # Batch lookup or filter-only - return list
            records = await execute_read(
                driver=driver,
                query=query,
                params=params,
                processor=process_multiple_records,
            )
            return [cls(**r["node"]) for r in records]
