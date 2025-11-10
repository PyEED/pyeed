from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class PyeedBase(BaseModel):
    """Base class for all nodes in the Database."""

    model_config = ConfigDict(frozen=False, validate_assignment=True, use_enum_values=True)

    custom: dict[str, Any] = Field(
        default_factory=dict, description="Arbitrary custom data as key-value pairs"
    )
