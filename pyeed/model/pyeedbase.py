import re
from collections.abc import Iterable as _Iter
from dataclasses import dataclass
from typing import Any, ClassVar, Dict, List, Optional, Tuple, Union

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
class EdgeMap:
    """Map of parent labels to relationship names."""

    rules: Dict[str, Union[str, Dict[str, str]]]

    def resolve(self, parent_label: str, field_name: Optional[str]) -> str:
        """Resolve the relationship name for a given parent + field.

        Args:
            parent_label: Name of the parent class/label.
            field_name: Field name of the parent, if applicable.

        Returns:
            The resolved relationship name.

        Raises:
            ValueError: If no mapping exists for the given parent/field.
        """
        spec = self.rules.get(parent_label)
        if spec is None:
            raise ValueError(
                f"No EdgeMap rule defined for parent '{parent_label}'. "
                f"Available parents: {list(self.rules.keys())}"
            )

        if isinstance(spec, str):
            return spec

        if field_name is None:
            raise ValueError(
                f"EdgeMap for parent '{parent_label}' requires a field name "
                f"(available: {list(spec.keys())})."
            )

        rel = spec.get(field_name)
        if rel is None:
            raise ValueError(
                f"No EdgeMap rule for field '{field_name}' under parent '{parent_label}'. "
                f"Available: {list(spec.keys())}"
            )

        return rel


class PyeedBase(BaseModel):
    """Base class for all nodes in the Database."""

    edge_map: ClassVar[Optional[EdgeMap]] = None

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
            parent_label = type(parent).__name__
            child_cls = type(child)
            emap = getattr(child_cls, "edge_map", None)

            if emap is None:
                # No mapping present: raise with a helpful message
                raise ValueError(
                    f"[Graphify] {child_cls.__name__} must define edge_map to connect "
                    f"from parent '{parent_label}' (field '{field_name}')."
                )

            # EdgeMap will raise a clear error if parent/field is not mapped
            rel_type = emap.resolve(parent_label=parent_label, field_name=field_name)

            sl, sk = parent_label, parent.get_unique_model_field()
            dl, dk = child_cls.__name__, child.get_unique_model_field()
            sv, dv = getattr(parent, sk), getattr(child, dk)

            edges.append(
                {
                    "type": rel_type,
                    "src": (sl, sk, sv),
                    "dst": (dl, dk, dv),
                }
            )

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
