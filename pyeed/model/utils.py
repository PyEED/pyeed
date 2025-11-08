from dataclasses import dataclass
from typing import get_args, get_origin

from .pyeedbase import LabelProperty, PyeedBase


@dataclass()
class UniqueSpec:
    label: str
    prop: str


@dataclass()
class BtreeSpec:
    label: str
    prop: str


@dataclass()
class VectorSpec:
    label: str
    prop: str


def collect_schema(
    models: list[type[PyeedBase]],
) -> tuple[list[UniqueSpec], list[BtreeSpec], list[VectorSpec]]:
    """Collect schema information from a list of models inheriting from BaseNode."""
    uniques: list[UniqueSpec] = []
    btrees: list[BtreeSpec] = []
    vectors: list[VectorSpec] = []

    for model in models:
        label = model.__name__
        for name, finfo in model.model_fields.items():
            ann = finfo.annotation
            if get_origin(ann) is not None:
                _, *meta = get_args(ann)
            else:
                meta = getattr(finfo, "metadata", [])

            for m in meta or ():
                if isinstance(m, LabelProperty):
                    if m.unique:
                        uniques.append(UniqueSpec(label, name))
                    if m.index:
                        btrees.append(BtreeSpec(label, name))
                    if m.vector_index:
                        vectors.append(VectorSpec(label, name))
    return uniques, btrees, vectors
