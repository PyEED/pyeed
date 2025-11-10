import re
from dataclasses import dataclass
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator


@dataclass(frozen=True)
class LabelProperty:
    """
    Propertys of a Neo4j label.
    Influences the Neo4j schema creation.
    """

    index: bool = False
    unique: bool = False


class PyeedBase(BaseModel):
    """Base class for all nodes in the Database."""

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
                    f"Nested dictionaries are not allowed in custom fields. Key '{key}' contains a dictionary value."
                )

        if conflicting_keys:
            raise ValueError(
                f"Custom field keys cannot conflict with existing attributes or be 'custom': {conflicting_keys}"
            )

        if invalid_keys:
            raise ValueError(
                f"Invalid custom field keys: {invalid_keys}. "
                f"Keys must start with a letter or underscore and contain only letters, digits, "
                f"or underscores."
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
        result = {k: v for k, v in flat.items() if _is_neo4j_prop_value(v)}
        return result

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
