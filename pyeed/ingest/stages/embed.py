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
from typing import TYPE_CHECKING, Any

from loguru import logger
from rich.progress import Progress, TaskID

from ...embed.esm2 import ESM2Embedder
from ...embed.types import distribute_embeddings_to_records, records_to_embedding_inputs
from ..core.pipeline import PipelineRecord
from ..core.protocol import SENTINEL, PipelineContext
from ..model.pyeedbase import PyeedBase

if TYPE_CHECKING:
    from ...db.milvus import VectorDB


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
        input_queue: asyncio.Queue[PipelineRecord[PyeedBase] | object] = next(
            iter(input_queues.values())
        )

        logger.info("Embedding stage starting")

        # Convert queue to record stream
        record_stream = self._queue_to_records(input_queue)

        # Process in chunks
        chunk_count = 0
        async for chunk in self._chunk_records(record_stream):
            # Check for total update on first item (pipeline is active by then)
            if progress is not None and task_id is not None:
                total = context.stats.get("total")
                if total is not None:
                    progress.update(task_id, total=total)
            chunk_count += 1
            logger.info(f"Processing chunk {chunk_count} with {len(chunk)} records")
            await self._process_chunk(chunk, output_queues, progress, task_id)

        # Send SENTINEL when done
        for queue in output_queues.values():
            await queue.put(SENTINEL)

        logger.info(f"Embedding stage finished. Processed {chunk_count} chunks")

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

    async def _records_to_batches(
        self,
        records: list[PipelineRecord[PyeedBase]],
    ) -> AsyncIterator[tuple[list[str], list[str]]]:
        """Convert records to length-sorted batches with natural backpressure.

        Args:
            records: List of pipeline records to convert

        Yields:
            (sequences, protein_ids) tuples ready for embedder
        """
        # Extract sequences and IDs from records
        sequences, protein_ids = records_to_embedding_inputs(records)

        # Create length-sorted batches for GPU efficiency
        batches = self.embedder.create_length_sorted_batches(
            sequences, protein_ids, self.batch_size
        )

        logger.debug(f"Created {len(batches)} batches from {len(records)} records")

        # Yield batches with async backpressure
        for batch in batches:
            await asyncio.sleep(0)  # Yield control to event loop for backpressure
            yield batch

    async def _process_chunk(
        self,
        records: list[PipelineRecord[PyeedBase]],
        output_queues: dict[str, asyncio.Queue[Any]],
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Process one chunk: check Milvus → batch → embed → convert → emit.

        Args:
            records: List of pipeline records to embed
            output_queues: Output queues to emit updated records
            progress: Progress instance for tracking
            task_id: Task ID for progress updates
        """
        if not records:
            logger.debug("Empty chunk, skipping")
            return

        # Check Milvus for existing embeddings (if configured)
        existing_ids: set[str] = set()
        if self.vector_db and self.collection_name:
            protein_ids = [record.data.id for record in records]
            try:
                existing_ids = await asyncio.to_thread(
                    self.vector_db._check_existing_ids,
                    protein_ids,
                )
                if existing_ids:
                    logger.info(
                        f"Found {len(existing_ids)} records with existing embeddings in Milvus, "
                        f"skipping embedding for those"
                    )
            except Exception as e:
                logger.warning(f"Failed to check existing IDs in Milvus: {e}")

        # Split: existing vs. new
        to_embed = [r for r in records if r.data.id not in existing_ids]
        already_embedded = [r for r in records if r.data.id in existing_ids]

        # Update progress for skipped records
        if progress is not None and task_id is not None and already_embedded:
            progress.update(task_id, advance=len(already_embedded))

        if already_embedded:
            for record in already_embedded:
                for queue in output_queues.values():
                    await queue.put(record)

        if not to_embed:
            logger.debug("No new records to embed after checking Milvus")
            return

        logger.debug(f"Embedding {len(to_embed)} new records")

        # Create async batch stream with natural backpressure
        batch_stream = self._records_to_batches(to_embed)

        # Stream through embedder - results arrive as they complete
        batch_count = 0
        async for embedding_batch in self.embedder.embed_stream(batch_stream):
            batch_count += 1
            logger.debug(
                f"Received embedding batch {batch_count} with "
                f"{len(embedding_batch.protein_ids)} proteins"
            )

            # Distribute embeddings to records
            distribute_embeddings_to_records(embedding_batch, to_embed)

            # Emit records that got embeddings (as they complete)
            batch_ids = set(embedding_batch.protein_ids)
            emitted = 0
            for record in to_embed:
                if record.data.id in batch_ids:
                    for queue in output_queues.values():
                        await queue.put(record)
                    emitted += 1

            logger.debug(f"Emitted {emitted} records to output queues")

            # Update progress
            if progress is not None and task_id is not None:
                progress.update(task_id, advance=len(embedding_batch.protein_ids))

        logger.debug(f"Chunk processing complete. Processed {batch_count} batches")
