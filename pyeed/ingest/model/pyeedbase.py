from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, TypeVar

from loguru import logger
from neo4j import AsyncDriver, AsyncTransaction
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
    def _check_indexed_fields(self):
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

        # cls._index_checked = True
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

        index_field_name = cls._get_indexed_field()

        cypher = f"""
        UNWIND $rows AS row
        MERGE (n:`{cls.__name__}` {{ {index_field_name}: row.{index_field_name} }})
        SET n += row
        """

        kwargs = session_kwargs or {}
        async with driver.session(**kwargs) as session:

            async def _tx(tx: AsyncTransaction, rows: list[dict[str, Any]]) -> None:
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
    def _get_indexed_field(cls) -> str:
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

    @classmethod
    def get_unique_field(cls) -> str:
        """Returns the name of the field marked with LabelProperty(index=True)."""
        return cls._get_indexed_field()


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


# ============================================================================
# TEST CODE - Remove before production
# ============================================================================
if __name__ == "__main__":
    from typing import Annotated

    from rich import print

    print("=" * 60)
    print("Testing BaseNode with LabelProperty(index=True)")
    print("=" * 60)

    # Test 1: Valid node with indexed field
    print("\n1. Testing valid node with indexed field:")
    try:

        class TestNode(BaseNode):
            iddd: Annotated[str, LabelProperty(index=True)] = Field(..., description="ID")
            name: str = Field(..., description="Name")

        node = TestNode(iddd="test-123", name="Test Node")
        print(f"   ✓ Created node: {node}")
        print(f"   ✓ unique_field: {node.get_unique_field()}")
        print(f"   ✓ unique_field (class): {TestNode.get_unique_field()}")
        print(f"   ✓ to_dict(): {node.to_dict()}")
    except Exception as e:
        print(f"   ✗ Error: {e}")

    # Test 2: Node without indexed field (should raise error)
    print("\n2. Testing node without indexed field (should fail):")
    try:

        class NoIndexNode(BaseNode):
            name: str = Field(..., description="Name")

        node = NoIndexNode(name="Test")
        print(f"   ✗ Should have failed! unique_field: {node.get_unique_field()}")
    except ValueError as e:
        print(f"   ✓ Correctly raised ValueError: {e}")

    # Test 3: Node with multiple indexed fields (should raise error)
    print("\n3. Testing node with multiple indexed fields (should fail):")
    try:

        class MultiIndexNode(BaseNode):
            id: Annotated[str, LabelProperty(index=True)] = Field(..., description="ID")
            code: Annotated[str, LabelProperty(index=True)] = Field(..., description="Code")

        node = MultiIndexNode(id="test", code="code")
        print(f"   ✗ Should have failed! Created: {node}")
    except ValueError as e:
        print(f"   ✓ Correctly raised ValueError: {e}")

    # Test 4: Custom field validation - dict values not allowed
    print("\n4. Testing custom field with dict value (should fail):")
    try:

        class TestNode2(BaseNode):
            id: Annotated[str, LabelProperty(index=True)] = Field(..., description="ID")

        node = TestNode2(id="test", custom={"nested": {"key": "value"}})
        print(f"   ✗ Should have failed! Created: {node}")
    except ValueError as e:
        print(f"   ✓ Correctly raised ValueError: {e}")

    # Test 5: Valid custom fields
    print("\n5. Testing valid custom fields:")
    try:

        class TestNode3(BaseNode):
            id: Annotated[str, LabelProperty(index=True)] = Field(..., description="ID")

        node = TestNode3(
            id="test", custom={"extra_field": "value", "number": 42, "tags": ["tag1", "tag2"]}
        )
        print(f"   ✓ Created node with custom fields: {node}")
        print(f"   ✓ to_dict(): {node.to_dict()}")
    except Exception as e:
        print(f"   ✗ Error: {e}")

    print("\n" + "=" * 60)
    print("Tests completed!")
    print("=" * 60)
