import numpy as np
import torch
from numpy.typing import NDArray

def mean_pooling(
    hidden_states: torch.Tensor, attention_mask: torch.Tensor | None = None
) -> torch.Tensor:
    """
    Apply mean pooling to hidden states.

    Args:
        hidden_states: Shape [batch, seq_len, hidden_dim]
        attention_mask: Shape [batch, seq_len], 1 for real tokens, 0 for padding

    Returns:
        Pooled embeddings [batch, hidden_dim]
    """
    if attention_mask is None:
        # Simple mean over sequence dimension
        return hidden_states.mean(dim=1)

    # Expand attention mask to match hidden_states dimensions
    # attention_mask: [batch, seq_len] -> [batch, seq_len, 1]
    expanded_mask = attention_mask.unsqueeze(-1).expand(hidden_states.size()).float()

    # Sum embeddings, weighted by mask
    sum_embeddings = torch.sum(hidden_states * expanded_mask, dim=1)

    # Sum mask values to get count of real tokens
    sum_mask = torch.clamp(expanded_mask.sum(dim=1), min=1e-9)

    # Average
    return sum_embeddings / sum_mask

def normalize_embeddings(embeddings: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    L2-normalize embeddings.

    Args:
        embeddings: Shape [batch, dim] or [dim]

    Returns:
        Normalized embeddings with same shape
    """
    if embeddings.ndim == 1:
        embeddings = embeddings.reshape(1, -1)
        squeeze = True
    else:
        squeeze = False

    norm = np.linalg.norm(embeddings, axis=1, keepdims=True)
    norm[norm == 0] = 1.0  # Avoid division by zero
    normalized = embeddings / norm

    if squeeze:
        return np.asarray(normalized[0], dtype=np.float64)
    return np.asarray(normalized, dtype=np.float64)
