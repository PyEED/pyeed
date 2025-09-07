import re
from dataclasses import dataclass
from typing import Any, ClassVar, Dict, List, Optional, Tuple

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
class ParentReference:
    """
    Contains information about the parent node of the current node.
    """

    parent_node_name: str
    rel_name: str


class PyeedBase(BaseModel):
    """Base class for all nodes in the Database."""

    PARENT_REF: ClassVar[Optional[ParentReference]] = None

    model_config = ConfigDict(
        frozen=False, validate_assignment=True, use_enum_values=True
    )

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
            current_class = (
                current_class.__bases__[0] if current_class.__bases__ else None
            )

        # Check for conflicts
        conflicting_keys = []
        invalid_keys = []

        # Regex pattern for valid Python variable names
        valid_var_pattern = re.compile(r"^[a-zA-Z_][a-zA-Z0-9_]*$")

        for key, value in v.items():
            # Check for conflicts with existing attributes
            if key in field_names:
                conflicting_keys.append(key)
            elif key == "custom":
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
                f"Custom field keys cannot conflict with existing attributes or be 'custom': {conflicting_keys}"
            )

        if invalid_keys:
            raise ValueError(
                f"Invalid custom field keys: {invalid_keys}. "
                f"Keys must start with a letter or underscore and contain only letters, digits, or underscores."
            )

        return v

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
            if (
                isinstance(field_info.metadata[0], LabelProperty)
                and field_info.metadata[0].unique
            ):
                return field_name

        raise ValueError(
            f"No unique field found. No field of {type(self)} is marked with NodeHint(unique=True)"
        )

    def graphify(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        nodes: List[Dict[str, Any]] = []
        edges: List[Dict[str, Any]] = []

        def emit_node(obj: "PyeedBase") -> Dict[str, Any]:
            lbl = obj.__class__.__name__
            k = obj.get_unique_model_field()
            v = getattr(obj, k)
            return {"label": lbl, "key": (k, v), "props": obj.to_dict()}

        def connect(parent: "PyeedBase", field_name: str, child: "PyeedBase") -> None:
            pref = getattr(type(child), "PARENT_REF", None) or getattr(
                type(child), "_parent_ref", None
            )
            # prefer PARENT_REF if it names this parent class
            if (
                pref
                and getattr(pref, "parent_node_name", None) == type(parent).__name__
            ):
                rel_type = getattr(pref, "rel_name", None) or getattr(
                    pref, "rel_type", None
                )
            else:
                rel_type = field_name.upper()

            sl, sk = type(parent).__name__, parent.get_unique_model_field()
            dl, dk = type(child).__name__, child.get_unique_model_field()
            sv, dv = getattr(parent, sk), getattr(child, dk)
            edges.append({"type": rel_type, "src": (sl, sk, sv), "dst": (dl, dk, dv)})

        from collections.abc import Iterable as _Iter

        seen: set[tuple[str, Any]] = set()

        def ensure_node(obj: "PyeedBase") -> None:
            """Ensure a node is added to the graph."""
            key = obj.get_unique_model_field()
            ident = (f"{type(obj).__name__}:{key}", getattr(obj, key))
            if ident not in seen:
                nodes.append(emit_node(obj))
                seen.add(ident)

        def walk(parent: "PyeedBase") -> None:
            """Walk the tree of nodes and edges."""
            for fname in parent.model_dump().keys():
                val = getattr(parent, fname, None)
                if isinstance(val, PyeedBase):
                    ensure_node(val)
                    connect(parent, fname, val)
                    walk(val)  # ← descend!
                elif isinstance(val, _Iter) and not isinstance(val, (str, bytes, dict)):
                    for item in val:
                        if isinstance(item, PyeedBase):
                            ensure_node(item)
                            connect(parent, fname, item)
                            walk(item)  # ← descend!

        ensure_node(self)
        walk(self)
        return nodes, edges


def _is_neo4j_primitive(x: object) -> bool:
    """Check if a value is a Neo4j primitive."""
    return isinstance(x, (str, int, float, bool)) or x is None


def _is_neo4j_prop_value(v: object) -> bool:
    """Check if a value is a Neo4j property value."""
    if _is_neo4j_primitive(v):
        return True
    if isinstance(v, (list, tuple)):
        return all(_is_neo4j_primitive(e) for e in v)
    return False
