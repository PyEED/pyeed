from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Generic, Self, TypeVar

from loguru import logger
from neo4j import AsyncDriver, AsyncTransaction
from neo4j._async.work import AsyncManagedTransaction
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

T = TypeVar("T", bound="BaseNode")


@dataclass(frozen=True)
class LabelProperty:
    """
    Influences the Neo4j schema creation.
    """

    index: bool = False


class BaseNode(BaseModel):
    """Base class for all nodes in the Database.

    To mark a field as indexed, use ``Annotated`` with ``LabelProperty(index=True)``
    as a type hint annotation. For example::

        from typing import Annotated
        from pydantic import Field

        class MyNode(BaseNode):
            id: Annotated[str, LabelProperty(index=True)] = Field(..., description="ID field")

    Note: Only one field per class can be marked with ``index=True``.
    """

    model_config = ConfigDict(frozen=False, validate_assignment=True, use_enum_values=True)

    custom: dict[str, Any] = Field(
        default_factory=dict, description="Arbitrary custom data as key-value pairs"
    )

    @model_validator(mode="after")
    def _check_indexed_fields(self) -> Self:
        cls = type(self)
        if getattr(cls, "_index_checked", False):
            return self

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

    @field_validator("custom")
    @classmethod
    def validate_keys(cls, v: dict[str, Any], info: Any) -> dict[str, Any]:
        """Validate that custom keys don't conflict with existing attributes."""
        if not v:
            return v

        # Get all field names from the current class and its parents
        field_names = set()
        current_class = info.data.get("__class__", cls)

        # Collect field names from current class and all parent classes
        while current_class and current_class != BaseModel:
            field_names.update(current_class.model_fields.keys())
            current_class = current_class.__bases__[0] if current_class.__bases__ else None

        # Check for conflicts
        conflicting_keys = []
        invalid_keys = []

        # Regex pattern for valid Python variable names
        valid_var_pattern = re.compile(r"^[a-zA-Z_][a-zA-Z0-9_]*$")

        for key, value in v.items():
            # Check for conflicts with existing attributes
            if key in field_names or key == "custom":
                conflicting_keys.append(key)

            # Check if key is a valid Python variable name
            if not valid_var_pattern.match(key):
                invalid_keys.append(key)

            # Check for nested dictionaries (not allowed)
            if isinstance(value, dict):
                raise ValueError(
                    f"Nested dictionaries are not allowed in custom fields. "
                    f"Key '{key}' contains a dictionary value."
                )

        if conflicting_keys:
            raise ValueError(
                f"Custom field keys cannot conflict with existing attributes or be 'custom': "
                f"{conflicting_keys}"
            )

        if invalid_keys:
            raise ValueError(
                f"Invalid custom field keys: {invalid_keys}. "
                f"Keys must start with a letter or underscore and contain only letters, digits, "
                f"or underscores."
            )

        return v

    @classmethod
    async def bulk_upsert(
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
        rows = [node.to_dict() for node in nodes]
        if not rows:
            return

        index_field_name = cls.get_index_field()

        cypher = f"""
        UNWIND $rows AS row
        MERGE (n:`{cls.__name__}` {{ {index_field_name}: row.{index_field_name} }})
        SET n += row
        """

        kwargs = session_kwargs or {}
        async with driver.session(**kwargs) as session:

            async def _tx(tx: AsyncManagedTransaction, rows: list[dict[str, Any]]) -> None:
                await tx.run(cypher, rows=rows)

            try:
                await session.execute_write(_tx, rows)
            except Exception:
                logger.exception("Failed to batch upsert nodes")

    async def upsert(
        self: T,
        driver: AsyncDriver,
        session_kwargs: dict[str, Any] | None = None,
    ) -> None:
        """Upsert this node instance.

        Args:
            driver: Neo4j async driver.
            session_kwargs: Optional session keyword arguments.
        """
        await type(self).bulk_upsert(
            driver,
            [self],
            session_kwargs=session_kwargs,
        )

    def to_dict(self) -> dict[str, Any]:
        """Convert the model to a Neo4j-safe dictionary.

        - flatten `custom`
        - keep only primitives or list-of-primitives
        """
        d = self.model_dump(exclude_none=True, exclude_unset=True)
        custom = d.pop("custom", {}) or {}
        flat = {**d, **custom}
        result = {k: v for k, v in flat.items() if _is_neo4j_prop_value(v)}
        return result

    @classmethod
    def get_index_field(cls) -> str:
        """Returns the name of the field marked with LabelProperty(index=True)."""
        for field_name, field_info in cls.model_fields.items():
            metadata = field_info.metadata
            for m in metadata:
                if isinstance(m, LabelProperty) and m.index:
                    return field_name

        raise ValueError(
            f"No field with LabelProperty(index=True) found in {cls.__name__}. "
            f"To mark a field as indexed, use: "
            f"Annotated[<type>, LabelProperty(index=True)] as the type hint annotation."
        )

def _is_primitive(x: object) -> bool:
    """Check if a value is a primitive type."""
    return isinstance(x, str | int | float | bool) or x is None


def _is_neo4j_prop_value(v: object) -> bool:
    """Check if a value is a Neo4j property value."""
    if _is_primitive(v):
        return True
    if isinstance(v, list | tuple):
        return all(_is_primitive(e) for e in v)
    return False



TStart = TypeVar("TStart", bound=BaseNode)
TEnd = TypeVar("TEnd", bound=BaseNode)


class RelationGroup(Generic[TStart, TEnd]):
    """One start node with multiple target nodes."""

    __slots__ = ("start", "targets")

    def __init__(self, start: TStart, targets: Iterable[TEnd]) -> None:
        self.start = start
        self.targets = list(targets)


RelationBatch = list[RelationGroup[TStart, TEnd]]


async def write_relations(
    driver: AsyncDriver,
    rel_type: str,
    batch: RelationBatch[TStart, TEnd],
    *,
    session_kwargs: dict[str, Any] | None = None,
    rel_props: dict[str, Any] | None = None,
) -> None:
    """Write (start)-[rel_type]->(target) relationships for a batch.

    All start nodes must be same type; all target nodes must be same type.
    """
    if not batch:
        return

    rows: list[dict[str, Any]] = []
    for group in batch:
        if not group.targets:
            continue
        from_id = getattr(group.start, group.start.get_index_field())
        from_label = group.start.__class__.__name__
        assert all(isinstance(t, (type(group.targets[0]))) for t in group.targets)
        for target in group.targets:
            to_id = getattr(target, target.get_index_field())
            to_label = target.__class__.__name__
            rows.append(
                {
                    "from_label": from_label,
                    "from_id": from_id,
                    "to_label": to_label,
                    "to_id": to_id,
                    "props": rel_props,
                }
            )

    if not rows:
        return

    #TODO: Check if the cypher syntax is valid

    cypher = f"""
    UNWIND $rows AS row
    MATCH (a:`{row.from_label}` {{ {row.from_id}: row.from_id }})
    MATCH (b:`{row.to_label}` {{ {row.to_id}: row.to_id }})
    MERGE (a)-[r:`{rel_type}`]->(b)
    SET r += row.props
    """

    session_kwargs = session_kwargs or {}
    async with driver.session(**session_kwargs) as session:

        async def _tx(tx: AsyncManagedTransaction, rows: list[dict[str, Any]]) -> None:
            await tx.run(cypher, rows=rows)

        try:
            await session.execute_write(_tx, rows)
        except Exception:
            logger.exception("Failed to write relations")
