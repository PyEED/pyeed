from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import torch

type TorchDTypeName = Literal["float32", "float16"]


@dataclass(slots=True, frozen=True)
class SequenceItem:
    id: str
    sequence: str


@dataclass(slots=True, frozen=True)
class TokenizedJob:
    """A job to be processed by the ESM2Processor."""

    batch: list[SequenceItem]
    input_ids: torch.Tensor  # (B, T)
    attention_mask: torch.Tensor  # (B, T)
    num_tokens: int


@dataclass(slots=True, frozen=True)
class EmbeddingRecord:
    """A record of an embedding."""

    id: str
    vector: np.ndarray
