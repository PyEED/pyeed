import re
from collections.abc import Iterable as _Iter
from dataclasses import dataclass
from typing import Any, ClassVar

from pydantic import BaseModel, ConfigDict, Field, field_validator


@dataclass(frozen=True)
class LabelProperty:
    """
    Propertys of a Neo4j label.
    Influences the Neo4j schema creation.
    """

    unique: bool = False
    index: bool = False
    vector_index: bool = False


@dataclass(frozen=True)
class Edge:
    parent_label: str
    rel_name: str
    field_name: str | None = None


class PyeedBase(BaseModel):
    """Base class for all nodes in the Database."""

    EDGES: ClassVar[tuple[Edge, ...]] = ()

    model_config = ConfigDict(frozen=False, validate_assignment=True, use_enum_values=True)

    custom: dict[str, Any] = Field(
        default_factory=dict, description="Arbitrary custom data as key-value pairs"
    )

    @field_validator("custom")
    @classmethod
    def validate_custom_keys(cls, v: dict[str, Any], info: Any) -> dict[str, Any]:
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
    def resolve_edge(cls, parent_label: str, field_name: str | None) -> str:
        exact = [
            e for e in cls.EDGES if e.parent_label == parent_label and e.field_name == field_name
        ]
        if len(exact) == 1:
            return exact[0].rel_name
        if len(exact) > 1:
            raise ValueError(
                f"{cls.__name__}: multiple edges match parent='{parent_label}', "
                f"field='{field_name}'."
            )
        generic = [e for e in cls.EDGES if e.parent_label == parent_label and e.field_name is None]
        if len(generic) == 1:
            return generic[0].rel_name
        if not generic:
            raise ValueError(
                f"{cls.__name__}: no edge for parent='{parent_label}' "
                f"(field='{field_name}'). Add to {cls.__name__}.EDGES."
            )
        raise ValueError(
            f"{cls.__name__}: ambiguous edges for parent='{parent_label}'. "
            f"Disambiguate by setting field_name."
        )

    def to_dict(self) -> dict[str, Any]:
        """Convert the model to a Neo4j-safe dictionary.

        - flatten `custom`
        - keep only primitives or list-of-primitives
        """
        d = self.model_dump(exclude_none=True, exclude_unset=True)
        custom = d.pop("custom", {}) or {}
        flat = {**d, **custom}
        return {k: v for k, v in flat.items() if _is_neo4j_prop_value(v)}

    def get_unique_model_field(self) -> str:
        """Returns the name of the field marked with NodeHint(unique=True)"""
        for field_name, field_info in type(self).model_fields.items():
            if not field_info.metadata:
                continue
            if isinstance(field_info.metadata[0], LabelProperty) and field_info.metadata[0].unique:
                return field_name

        raise ValueError(
            f"No unique field found. No field of {type(self)} is marked with NodeHint(unique=True)"
        )

    def graphify(self) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
        nodes: list[dict[str, Any]] = []
        rels: list[dict[str, Any]] = []
        seen: set[tuple[str, Any]] = set()

        def emit_node(obj: "PyeedBase") -> dict[str, Any]:
            lbl = obj.__class__.__name__
            k = obj.get_unique_model_field()
            v = getattr(obj, k)
            return {"label": lbl, "key": (k, v), "props": obj.to_dict()}

        def ensure_node(obj: "PyeedBase") -> None:
            key = obj.get_unique_model_field()
            ident = (f"{type(obj).__name__}:{key}", getattr(obj, key))
            if ident not in seen:
                nodes.append(emit_node(obj))
                seen.add(ident)

        def connect(parent: "PyeedBase", field_name: str, child: "PyeedBase") -> None:
            parent_label = type(parent).__name__
            child_cls = type(child)

            if not hasattr(child_cls, "EDGES"):
                raise ValueError(
                    f"{child_cls.__name__} must define EDGES "
                    f"to connect from parent '{parent_label}' via field '{field_name}'."
                )

            rel_type = child_cls.resolve_edge(parent_label=parent_label, field_name=field_name)

            sl, sk = parent_label, parent.get_unique_model_field()
            dl, dk = child_cls.__name__, child.get_unique_model_field()
            sv, dv = getattr(parent, sk), getattr(child, dk)

            rels.append({"type": rel_type, "src": (sl, sk, sv), "dst": (dl, dk, dv)})

        def walk(parent: "PyeedBase") -> None:
            for fname in parent.model_dump():
                val = getattr(parent, fname, None)
                if isinstance(val, PyeedBase):
                    ensure_node(val)
                    connect(parent, fname, val)
                    walk(val)
                elif isinstance(val, _Iter) and not isinstance(val, str | bytes | dict):
                    for item in val:
                        if isinstance(item, PyeedBase):
                            ensure_node(item)
                            connect(parent, fname, item)
                            walk(item)

        ensure_node(self)
        walk(self)
        return nodes, rels


def _is_neo4j_primitive(x: object) -> bool:
    """Check if a value is a Neo4j primitive."""
    return isinstance(x, str | int | float | bool) or x is None


def _is_neo4j_prop_value(v: object) -> bool:
    """Check if a value is a Neo4j property value."""
    if _is_neo4j_primitive(v):
        return True
    if isinstance(v, list | tuple):
        return all(_is_neo4j_primitive(e) for e in v)
    return False
