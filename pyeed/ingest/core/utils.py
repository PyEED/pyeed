"""Reusable utilities for pipeline stages."""

from __future__ import annotations

import asyncio
from collections import defaultdict

from pyeed.ingest.core.pipeline import IngestItem
from pyeed.ingest.model.pyeedbase import BaseNode

# async def queue_to_async_iter(
#     queue: asyncio.Queue[T | object],
# ) -> AsyncIterator[T]:
#     """Convert a queue to an async iterator that stops at SENTINEL.

#     Args:
#         queue: Queue to consume from

#     Yields:
#         Items from queue until SENTINEL is received

#     Example:
#         >>> async for item in queue_to_async_iter(my_queue):
#         ...     process(item)
#     """
#     while True:
#         item = await queue.get()
#         if item is SENTINEL:
#             break
#         yield item  # type: ignore


# async def batch_accumulator(
#     source: AsyncIterator[T],
#     batch_size: int,
# ) -> AsyncIterator[list[T]]:
#     """Accumulate items from source into batches.

#     Args:
#         source: Async iterator to batch
#         batch_size: Maximum items per batch

#     Yields:
#         Lists of items (last batch may be smaller)

#     Example:
#         >>> async for batch in batch_accumulator(items, batch_size=32):
#         ...     process_batch(batch)
#     """
#     batch: list[T] = []
#     async for item in source:
#         batch.append(item)
#         if len(batch) >= batch_size:
#             yield batch
#             batch = []

#     if batch:  # Final partial batch
#         yield batch


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


def group_by_relation(
    items: list[IngestItem[BaseNode]],
    relation_type: str,
) -> dict[str, list[BaseNode]]:
    """Group items by relation type.

    Args:
        items: List of ingest items.
        relation_type: Relation type to group by.

    Returns:
        Dictionary of relation type to list of nodes.
    """
    rels = defaultdict[str, list[BaseNode]](list)
    for item in items:
        if not item.relations:
            continue

        for rel_type, rel_values in item.relations.items():
            if rel_type != relation_type:
                continue

            for rel_value in rel_values:
                rels[rel_value].append(item.node)

    return rels
