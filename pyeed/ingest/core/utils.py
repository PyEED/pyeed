"""Reusable utilities for pipeline stages."""

from __future__ import annotations

import asyncio
from collections.abc import AsyncIterator
from typing import TypeVar

from .protocol import SENTINEL

type T = TypeVar("T")


async def queue_to_async_iter(
    queue: asyncio.Queue[T | object],
) -> AsyncIterator[T]:
    """Convert a queue to an async iterator that stops at SENTINEL.

    Args:
        queue: Queue to consume from

    Yields:
        Items from queue until SENTINEL is received

    Example:
        >>> async for item in queue_to_async_iter(my_queue):
        ...     process(item)
    """
    while True:
        item = await queue.get()
        if item is SENTINEL:
            break
        yield item  # type: ignore


async def batch_accumulator(
    source: AsyncIterator[T],
    batch_size: int,
) -> AsyncIterator[list[T]]:
    """Accumulate items from source into batches.

    Args:
        source: Async iterator to batch
        batch_size: Maximum items per batch

    Yields:
        Lists of items (last batch may be smaller)

    Example:
        >>> async for batch in batch_accumulator(items, batch_size=32):
        ...     process_batch(batch)
    """
    batch: list[T] = []
    async for item in source:
        batch.append(item)
        if len(batch) >= batch_size:
            yield batch
            batch = []

    if batch:  # Final partial batch
        yield batch


async def sort_by_length_batch(
    batch: list[tuple[str, str]],
) -> list[tuple[str, str]]:
    """Sort (id, sequence) pairs by sequence length (descending).

    Runs in thread pool to avoid blocking event loop.
    """
    return await asyncio.to_thread(
        sorted,
        batch,
        key=lambda x: len(x[1]),
        reverse=True,
    )
