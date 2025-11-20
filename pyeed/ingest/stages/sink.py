"""Neo4j and Milvus sink stages for batch processing PipelineRecords."""

from __future__ import annotations

import asyncio
from collections import defaultdict
from typing import Any

from loguru import logger
from rich.progress import Progress, TaskID

from pyeed.db.milvus import VectorDB
from pyeed.db.neo4j import GraphDB
from pyeed.db.queries import upsert_pipeline_records
from pyeed.embed.types import EmbeddingBatch
from pyeed.ingest.core.pipeline import PipelineContext, PipelineRecord
from pyeed.ingest.core.protocol import SENTINEL
from pyeed.ingest.model.pyeedbase import PyeedBase


class Neo4jUpsertStage:
    """Batches PipelineRecord objects from queue and upserts to Neo4j.

    Processes PipelineRecords with their children and relationships in batches
    to minimize database round trips.
    """

    def __init__(self, db: GraphDB, batch_size: int):
        """Initialize Neo4j sink.

        Args:
            db: GraphDB instance
            batch_size: Number of records to accumulate before upserting
        """
        self.db = db
        self.batch_size = batch_size

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Consume PipelineRecords, batch, upsert to Neo4j, and forward to embedding."""
        input_queue: asyncio.Queue[PipelineRecord[PyeedBase]] = next(iter(input_queues.values()))
        output_queue = next(iter(output_queues.values())) if output_queues else None

        item_dict: dict[str, set[str]] = defaultdict(set)
        batch: list[PipelineRecord[PyeedBase]] = []

        while True:
            item = await input_queue.get()

            if item is SENTINEL:
                # Flush remaining batch
                if batch:
                    async with self.db.async_driver.session() as session:
                        await upsert_pipeline_records(session, batch)
                    # Track and forward
                    for record in batch:
                        unique_field = record.data.get_unique_model_field()
                        unique_value = str(getattr(record.data, unique_field))
                        context.added_nodes[type(record.data).__name__].add(unique_value)
                        if output_queue:
                            await output_queue.put(record)
                    if progress is not None and task_id is not None:
                        progress.advance(task_id, len(batch))

                # Forward SENTINEL
                if output_queue:
                    await output_queue.put(SENTINEL)
                break

            unique_field = item.data.get_unique_model_field()
            unique_value = str(getattr(item.data, unique_field))
            item_dict[type(item.data).__name__].add(unique_value)

            batch.append(item)

            if len(batch) >= self.batch_size:
                async with self.db.async_driver.session() as session:
                    await upsert_pipeline_records(session, batch)
                # Track and forward
                for record in batch:
                    unique_field = record.data.get_unique_model_field()
                    unique_value = str(getattr(record.data, unique_field))
                    context.added_nodes[type(record.data).__name__].add(unique_value)
                    if output_queue:
                        await output_queue.put(record)
                if progress is not None and task_id is not None:
                    progress.advance(task_id, len(batch))
                batch = []


class MilvusUpsertStage:
    """Batches PipelineRecord embeddings and upserts to Milvus.

    Only processes records with embeddings. Sequences are NOT stored
    in Milvus (already in Neo4j).
    """

    def __init__(
        self,
        vector_db: VectorDB,
        collection_name: str,
        batch_size: int = 100,
    ):
        """Initialize Milvus sink.

        Args:
            vector_db: VectorDB instance
            collection_name: Milvus collection name
            batch_size: Records to accumulate before writing
        """
        self.vector_db = vector_db
        self.collection_name = collection_name
        self.batch_size = batch_size

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Consume PipelineRecords with embeddings, batch, upsert to Milvus."""
        input_queue = next(iter(input_queues.values()))
        batch: list[PipelineRecord[PyeedBase]] = []

        while True:
            item = await input_queue.get()

            if item is SENTINEL:
                # Flush remaining batch
                if batch:
                    embedding_batch = self._convert_to_embedding_batch(batch)
                    inserted = await self.vector_db.insert_async(
                        self.collection_name,
                        embedding_batch,
                        include_sequence=False,  # Don't store sequences
                    )
                    if progress is not None and task_id is not None:
                        progress.advance(task_id, inserted)
                break

            # Skip records without embeddings
            if not item.embeddings:
                logger.debug(f"Skipping record without embeddings: {item.data.id}")
                if progress is not None and task_id is not None:
                    progress.advance(task_id, 1)
                continue

            batch.append(item)

            if len(batch) >= self.batch_size:
                embedding_batch = self._convert_to_embedding_batch(batch)
                inserted = await self.vector_db.insert_async(
                    self.collection_name,
                    embedding_batch,
                    include_sequence=False,  # Don't store sequences
                )
                if progress is not None and task_id is not None:
                    progress.advance(task_id, inserted)
                batch = []

    def _convert_to_embedding_batch(
        self,
        records: list[PipelineRecord[PyeedBase]],
    ) -> EmbeddingBatch:
        """Convert PipelineRecords to EmbeddingBatch.

        Args:
            records: List of PipelineRecord with embeddings

        Returns:
            EmbeddingBatch with combined embeddings
        """
        batch = EmbeddingBatch()

        for record in records:
            # Extract protein ID and sequence
            batch.protein_ids.append(record.data.id)
            batch.sequences.append(record.data.sequence)

            # Group embeddings by pooling method
            for pooling_name, embedding in record.embeddings.items():
                if pooling_name not in batch.embeddings:
                    batch.embeddings[pooling_name] = []
                batch.embeddings[pooling_name].append(embedding)

        return batch
