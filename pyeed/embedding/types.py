"""Type definitions for embeddings."""

from typing import Literal

import numpy as np
import torch
from numpy.typing import NDArray

# Type literals for model and return dtypes
type ModelDType = Literal["float16", "bfloat16", "float32"]
type ReturnDType = Literal["float16", "bfloat16", "float32"]

# Dtype mappings
TF_DTYPE_MAP: dict[ModelDType, torch.dtype] = {
    "float16": torch.float16,
    "float32": torch.float32,
}

NP_DTYPE_MAP: dict[ReturnDType, type[np.generic]] = {
    "float16": np.float16,
    "float32": np.float32,
}

# Type alias for embedding arrays (float types only)
type EmbeddingArray = NDArray[np.float32 | np.float16]
