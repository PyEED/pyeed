"""Embedding sink protocols and implementations for storing embeddings."""

from __future__ import annotations

from typing import Any, Protocol

__all__ = [
    "EmbeddingSink",
]


# Type alias for embedding results (imported dynamically to avoid circular imports)
# EmbeddingResult = dict[str, list[EmbeddingArray] | list[str]]
# where EmbeddingArray = NDArray[np.float32 | np.float16]


class EmbeddingSink(Protocol):
    """Protocol for embedding storage backends.

    Implementations must handle streaming embedding batches and provide
    lifecycle management (connect, disconnect, flush).
    """

    async def write_batch(
        self,
        result: dict[str, Any],
        metadata: dict[str, Any] | None = None,
    ) -> int:
        """Write a batch of embeddings to the sink.

        Args:
            result: Dict mapping pooling method names to embedding arrays.
                   Must contain "ids" key with sequence identifiers.
            metadata: Optional metadata (model_name, layer_index, etc.)

        Returns:
            Number of embeddings written.
        """
        ...

    async def connect(self) -> None:
        """Initialize connection to the sink."""
        ...

    async def disconnect(self) -> None:
        """Close connection and flush any buffered data."""
        ...

    async def flush(self) -> None:
        """Flush any buffered writes."""
        ...
