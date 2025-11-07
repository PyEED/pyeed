"""Type definitions for embeddings."""

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import torch
from numpy.typing import NDArray

# Type literals for model and return dtypes
type ModelDType = Literal["float16", "bfloat16", "float32"]
type ReturnDType = Literal["float16", "bfloat16", "float32"]


type Vector = NDArray[np.float32] | NDArray[np.float16]


# Dtype mappings
TF_DTYPE_MAP: dict[ModelDType, torch.dtype] = {
    "float16": torch.float16,
    "float32": torch.float32,
}

NP_DTYPE_MAP: dict[ReturnDType, type[np.generic]] = {
    "float16": np.float16,
    "float32": np.float32,
}


@dataclass(slots=True)
class EmbeddingBatch:
    protein_ids: list[str] = field(default_factory=list)
    sequences: list[str] = field(default_factory=list)
    embeddings: dict[str, list[Vector]] = field(default_factory=dict)


@dataclass(slots=True)
class EmbeddingRecord:
    protein_id: str
    sequence: str
    embedding: dict[str, Vector]


def flatten_embedding_record(
    record: EmbeddingRecord,
) -> dict[str, str | Vector]:
    """Flatten an embedding record into a dictionary."""

    d: dict[str, str | Vector] = {"protein_id": record.protein_id, "sequence": record.sequence}
    d.update(record.embedding)
    return d
