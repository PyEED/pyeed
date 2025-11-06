from collections.abc import Sequence
from typing import Protocol, runtime_checkable

import numpy as np
import torch


@runtime_checkable
class PoolingFn(Protocol):
    def __call__(
        self,
        hidden_states: torch.Tensor,  # [B, L, D]
        attention_mask: torch.Tensor | None,  # [B, L] or None
    ) -> torch.Tensor:  # [B, D]
        ...


PoolingLike = PoolingFn | None | Sequence[PoolingFn | None]


def l2_normalize(
    x: torch.Tensor,
    dim: int = -1,
    eps: float = 1e-12,
) -> torch.Tensor:
    orig_dtype = x.dtype

    # accumulate norm in fp32 for stability
    xt32 = x if orig_dtype == torch.float32 else x.to(torch.float32)
    norm = torch.linalg.vector_norm(xt32, ord=2, dim=dim, keepdim=True).clamp_min(eps)
    y32 = xt32 / norm
    y = y32.to(orig_dtype) if orig_dtype != torch.float32 else y32

    return y


def normalize_cast_renorm(
    x: np.ndarray | torch.Tensor,
    cast_dtype: torch.dtype,
    dim: int = -1,
    eps: float = 1e-12,
) -> np.ndarray:
    y = l2_normalize(x, dim=dim, eps=eps)  # fp32 accumulation
    y = y.to(cast_dtype)
    y = l2_normalize(y, dim=dim, eps=eps)  # fix norm drift after cast

    return y


def mean_pooling(
    hidden_states: torch.Tensor,  # [B, L, D]
    attention_mask: torch.Tensor | None = None,  # [B, L] or None
) -> torch.Tensor:  # [B, D]
    """
    Mean pooling over sequence dimension.

    Args:
        hidden_states: [batch, seq_len, dim]
        attention_mask: [batch, seq_len], 1=real token, 0=padding

    Returns:
        [batch, dim] pooled embeddings
    """
    if attention_mask is None:
        pooled = hidden_states.mean(dim=1)
    else:
        # Expand [B, L] -> [B, L, 1] for broadcasting, no dtype conversion
        mask = attention_mask.unsqueeze(-1)  # [B, L, 1]

        # Multiply and sum; dtype determined by hidden_states
        summed = (hidden_states * mask).sum(dim=1)  # [B, D]
        counts = mask.sum(dim=1)  # [B, 1]
        pooled = summed / counts  # [B, D]

    return pooled


def max_pooling(
    hidden_states: torch.Tensor,  # [B, L, D]
    attention_mask: torch.Tensor | None = None,  # [B, L] or None
) -> torch.Tensor:  # [B, D]
    """
    Max pooling over sequence dimension.

    Args:
        hidden_states: [batch, seq_len, dim]
        attention_mask: [batch, seq_len], 1=real token, 0=padding

    Returns:
        [batch, dim] pooled embeddings
    """
    if attention_mask is None:
        pooled = hidden_states.max(dim=1).values
    else:
        # Expand [B, L] -> [B, L, 1] for broadcasting
        mask = attention_mask.unsqueeze(-1)  # [B, L, 1]

        # Fill padded positions with -inf so they don't affect max
        neg_inf = torch.finfo(hidden_states.dtype).min
        masked = torch.where(mask.bool(), hidden_states, neg_inf)

        pooled = masked.max(dim=1).values  # [B, D]

    return pooled


def min_pooling(
    hidden_states: torch.Tensor,  # [B, L, D]
    attention_mask: torch.Tensor | None = None,  # [B, L] or None
) -> torch.Tensor:  # [B, D]
    """
    Min pooling over sequence dimension.

    Args:
        hidden_states: [batch, seq_len, dim]
        attention_mask: [batch, seq_len], 1=real token, 0=padding

    Returns:
        [batch, dim] pooled embeddings
    """
    if attention_mask is None:
        pooled = hidden_states.min(dim=1).values
    else:
        # Expand [B, L] -> [B, L, 1] for broadcasting
        mask = attention_mask.unsqueeze(-1)  # [B, L, 1]

        # Fill padded positions with +inf so they don't affect min
        pos_inf = torch.finfo(hidden_states.dtype).max
        masked = torch.where(mask.bool(), hidden_states, pos_inf)

        pooled = masked.min(dim=1).values  # [B, D]

    return pooled
