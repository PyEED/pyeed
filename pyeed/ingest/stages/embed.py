"""Embedding stage for pipeline processing.

Streams records from queue, batches them, and processes through ESM2 embedder
using a clean streaming architecture with natural backpressure.

Architecture:
    Queue → Iterator → Chunker → Batcher → Embedder → Distributor → Queue

Backpressure:
    When embedder GPU tasks saturate (max concurrent limit reached), batches
    accumulate at the embedder. This blocks batch streaming, which fills the
    pipeline queue, which stops the reader from ingesting more sequences.

The embedder uses task-based parallelism (no internal queues) with a sliding
window of concurrent GPU tasks for efficient parallel processing while maintaining
natural backpressure to the pipeline.
"""

from __future__ import annotations

import asyncio
from collections.abc import AsyncIterator
from typing import TYPE_CHECKING, Any, TypeAlias

from loguru import logger
from rich.progress import Progress, TaskID

from ...embed.esm2 import ESM2Embedder
from ...embed.types import (
    EmbeddingBatch,
    distribute_embeddings_to_records,
    records_to_embedding_inputs,
)
from ..core.pipeline import PipelineRecord
from ..core.protocol import SENTINEL, PipelineContext

if TYPE_CHECKING:
    from ...db.milvus import VectorDB

# Type alias for clarity
Batch: TypeAlias = tuple[list[str], list[str]]


class EmbeddingStage:
    """Stream records, batch, embed, and emit results.

    Processes records in chunks with length-sorted batching for GPU efficiency.
    Uses async iterators throughout for composability and natural backpressure.
    """

    def __init__(
        self,
        embedder: ESM2Embedder,
        chunk_size: int = 1000,
        batch_size: int = 32,
        vector_db: VectorDB | None = None,
        collection_name: str | None = None,
    ):
        """Initialize embedding stage.

        Args:
            embedder: ESM2Embedder instance (must be initialized before pipeline run)
            chunk_size: Number of records to accumulate before embedding (default: 1000)
            batch_size: GPU batch size for embedding (default: 32)
            vector_db: Optional VectorDB to check for existing embeddings
            collection_name: Collection name to check
        """
        self.embedder = embedder
        self.chunk_size = chunk_size
        self.batch_size = batch_size
        self.vector_db = vector_db
        self.collection_name = collection_name

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Main processing loop.

        Converts input queue to record stream, processes in chunks, and streams results.

        Args:
            input_queues: Dictionary of input queues (expects single queue)
            output_queues: Dictionary of output queues
            context: Shared pipeline context
            progress: Progress instance for tracking
            task_id: Task ID for progress updates
        """
        input_queue: asyncio.Queue[PipelineRecord[BaseNode] | object] = next(
            iter(input_queues.values())
        )

        logger.info("Embedding stage starting")

        # Bounded queue to hold length-sorted batches ready for GPU consumption.
        # Size ties to number of devices to preserve backpressure upstream.
        prepared_batch_queue: asyncio.Queue[Batch | object] = asyncio.Queue(
            maxsize=max(len(self.embedder.devices) * 4, 4)
        )
        record_map: dict[str, PipelineRecord[PyeedBase]] = {}
        record_map_lock = asyncio.Lock()

        async def batch_producer() -> None:
            """Read from input queue, dedupe via Milvus, length-sort, enqueue batches."""
            chunk_count = 0
            async for chunk in self._chunk_records(self._queue_to_records(input_queue)):
                # Check for total update on first item (pipeline is active by then)
                if progress is not None and task_id is not None:
                    total = context.stats.get("total")
                    if total is not None:
                        progress.update(task_id, total=total)

                chunk_count += 1
                logger.info(f"Preparing chunk {chunk_count} with {len(chunk)} records")

                to_embed, already_embedded = await self._filter_existing(chunk)

                # Update progress and forward already embedded records immediately
                if progress is not None and task_id is not None and already_embedded:
                    progress.update(task_id, advance=len(already_embedded))

                if already_embedded:
                    for record in already_embedded:
                        for queue in output_queues.values():
                            await queue.put(record)

                if not to_embed:
                    logger.debug("Chunk had no new records after dedupe; skipping embedding")
                    continue

                sequences, protein_ids = records_to_embedding_inputs(to_embed)

                # Track records so consumer can attach embeddings
                async with record_map_lock:
                    for record in to_embed:
                        record_map[record.data.id] = record

                batches = self.embedder.create_length_sorted_batches(
                    sequences,
                    protein_ids,
                    self.batch_size,
                )
                logger.debug(f"Enqueuing {len(batches)} batches from chunk {chunk_count}")

                for batch in batches:
                    await prepared_batch_queue.put(batch)

            # Signal completion downstream
            await prepared_batch_queue.put(SENTINEL)
            logger.info(f"Batch producer finished after {chunk_count} chunks")

        async def batch_consumer() -> None:
            """Drain prepared batches, run embedding, distribute, emit."""
            batch_id = 0

            async def prepared_batch_stream() -> AsyncIterator[Batch]:
                while True:
                    item = await prepared_batch_queue.get()
                    if item is SENTINEL:
                        # Sentinel consumed; do not re-enqueue to avoid double-stop
                        return
                    yield item

            async for embedding_batch in self.embedder.embed_stream(prepared_batch_stream()):
                batch_id += 1
                await self._emit_embedded_records(
                    embedding_batch,
                    record_map,
                    record_map_lock,
                    output_queues,
                    progress,
                    task_id,
                )

            logger.info(f"Batch consumer finished after processing {batch_id} batches")

        # Launch producer/consumer concurrently to overlap preparation with GPU work
        await asyncio.gather(batch_producer(), batch_consumer())

        # Send SENTINEL when done to downstream stages
        for queue in output_queues.values():
            await queue.put(SENTINEL)

        logger.info("Embedding stage finished")

    async def _queue_to_records(
        self, queue: asyncio.Queue[PipelineRecord[PyeedBase] | object]
    ) -> AsyncIterator[PipelineRecord[PyeedBase]]:
        """Convert queue to async iterator, stopping at SENTINEL.

        Args:
            queue: Input queue containing PipelineRecord objects

        Yields:
            PipelineRecord objects from the queue
        """
        while True:
            item = await queue.get()
            if item is SENTINEL:
                logger.debug("Received SENTINEL in record stream")
                break
            yield item  # type: ignore

    async def _chunk_records(
        self,
        stream: AsyncIterator[PipelineRecord[PyeedBase]],
    ) -> AsyncIterator[list[PipelineRecord[PyeedBase]]]:
        """Accumulate records from stream into chunks.

        Args:
            stream: Async iterator of PipelineRecord objects

        Yields:
            Lists of records (chunks) of size self.chunk_size
        """
        chunk: list[PipelineRecord[PyeedBase]] = []
        async for record in stream:
            chunk.append(record)
            if len(chunk) >= self.chunk_size:
                yield chunk
                chunk = []

        # Yield remaining records
        if chunk:
            yield chunk

    async def _filter_existing(
        self, records: list[PipelineRecord[PyeedBase]]
    ) -> tuple[list[PipelineRecord[PyeedBase]], list[PipelineRecord[PyeedBase]]]:
        """Split records into new vs existing based on Milvus presence."""
        if not records or not (self.vector_db and self.collection_name):
            return records, []

        protein_ids = [record.data.id for record in records]
        try:
            existing_ids = await asyncio.to_thread(
                self.vector_db._check_existing_ids,
                protein_ids,
            )
            if existing_ids:
                logger.info(
                    "Skipping already embedded records",
                    extra={"existing": len(existing_ids)},
                )
        except Exception as e:
            logger.warning(f"Failed to check existing IDs in Milvus: {e}")
            existing_ids = set()

        to_embed = [r for r in records if r.data.id not in existing_ids]
        already_embedded = [r for r in records if r.data.id in existing_ids]
        return to_embed, already_embedded

    async def _emit_embedded_records(
        self,
        embedding_batch: EmbeddingBatch,
        record_map: dict[str, PipelineRecord[PyeedBase]],
        record_map_lock: asyncio.Lock,
        output_queues: dict[str, asyncio.Queue[Any]],
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Attach embeddings to records and emit downstream."""
        async with record_map_lock:
            try:
                records = [record_map[pid] for pid in embedding_batch.protein_ids]
            except KeyError as err:
                missing = err.args[0]
                raise ValueError(f"EmbeddingBatch contains unknown protein_id {missing}") from err

            distribute_embeddings_to_records(embedding_batch, records)

            for record in records:
                record_map.pop(record.data.id, None)

        for record in records:
            for queue in output_queues.values():
                await queue.put(record)

        if progress is not None and task_id is not None:
            progress.update(task_id, advance=len(embedding_batch.protein_ids))
