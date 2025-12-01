from collections.abc import Iterable
from dataclasses import dataclass

from .pyeedbase import BaseNode, LabelProperty


@dataclass(frozen=True, slots=True)
class IndexSpec:
    """Index specification."""

    label: str
    prop: str


def collect_schema(models: Iterable[type[BaseNode]]) -> list[IndexSpec]:
    """Collect the schema of the models and return the indices."""
    indices: list[IndexSpec] = []
    for model in models:
        label = model.__name__
        for name, finfo in model.model_fields.items():
            for meta in getattr(finfo, "metadata", ()):
                if not isinstance(meta, LabelProperty):
                    continue
                if meta.index:
                    indices.append(IndexSpec(label, name))

    return indices
