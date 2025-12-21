"""Neo4j and Milvus sink stages for batch processing PipelineRecords.

Note: MilvusUpsertStage is deprecated. Use pyeed.db.milvus.bulk_insert() with
the new async Milvus client for embedding storage.
"""

from __future__ import annotations

import asyncio
import warnings
from collections import defaultdict
from typing import TYPE_CHECKING, Any

from loguru import logger
from neo4j import AsyncDriver
from rich.progress import Progress, TaskID

from pyeed.db.queries import query_existing_nodes_by_ids, upsert_pipeline_records
from pyeed.ingest.core.pipeline import PipelineContext, PipelineRecord
from pyeed.ingest.core.protocol import SENTINEL
from pyeed.ingest.model.pyeedbase import BaseNode

if TYPE_CHECKING:
    from typing import Protocol

    class VectorDBProtocol(Protocol):
        def _check_existing_ids(self, ids: list[str]) -> set[str]: ...
        async def insert_async(self, collection: str, batch: Any, **kw: object) -> None: ...


class Neo4jUpsertStage:
    """Batches PipelineRecord objects from queue and upserts to Neo4j.

    Processes PipelineRecords with their children and relationships in batches
    to minimize database round trips.
    """

    def __init__(self, driver: AsyncDriver, batch_size: int):
        """Initialize Neo4j sink.

        Args:
            driver: Neo4j async driver
            batch_size: Number of records to accumulate before upserting
        """
        self.driver = driver
        self.batch_size = batch_size

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Consume PipelineRecords, batch, upsert to Neo4j, and forward to all output queues."""
        input_queue: asyncio.Queue[PipelineRecord[BaseNode]] = next(iter(input_queues.values()))

        item_dict: dict[str, set[str]] = defaultdict(set)
        batch: list[PipelineRecord[BaseNode]] = []

        while True:
            item = await input_queue.get()

            # Check for total update on first item (pipeline is active by then)
            if progress is not None and task_id is not None:
                total = context.stats.get("total")
                if total is not None:
                    progress.update(task_id, total=total)

            if item is SENTINEL:
                # Flush remaining batch
                if batch:
                    await self._process_batch(batch, output_queues, context, progress, task_id)
                # Forward SENTINEL to ALL output queues
                for output_queue in output_queues.values():
                    await output_queue.put(SENTINEL)
                break

            unique_field = item.data.get_unique_model_field()
            unique_value = str(getattr(item.data, unique_field))
            item_dict[type(item.data).__name__].add(unique_value)

            batch.append(item)

            if len(batch) >= self.batch_size:
                await self._process_batch(batch, output_queues, context, progress, task_id)
                batch = []

    async def _process_batch(
        self,
        batch: list[PipelineRecord[BaseNode]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Process a batch: check existing nodes, filter, upsert new ones, forward only new ones."""
        if not batch:
            return

        # Group records by label for batch existence checking
        by_label: dict[str, list[PipelineRecord[BaseNode]]] = defaultdict(list)
        for record in batch:
            label = type(record.data).__name__
            by_label[label].append(record)

        # Check which nodes already exist, grouped by label
        existing_ids: set[tuple[str, str]] = set()  # (label, unique_value)
        async with self.driver.session() as session:
            for label, label_records in by_label.items():
                if not label_records:
                    continue

                unique_field = label_records[0].data.get_unique_model_field()
                unique_values = [
                    str(getattr(record.data, unique_field)) for record in label_records
                ]

                existing_unique_values = await query_existing_nodes_by_ids(
                    session, label, unique_field, unique_values
                )

                for unique_value in existing_unique_values:
                    existing_ids.add((label, unique_value))

        # Filter batch: keep only records where parent doesn't exist
        new_batch: list[PipelineRecord[BaseNode]] = []
        for record in batch:
            label = type(record.data).__name__
            unique_field = record.data.get_unique_model_field()
            unique_value = str(getattr(record.data, unique_field))

            if (label, unique_value) not in existing_ids:
                new_batch.append(record)

        # Upsert only new records
        if new_batch:
            async with self.driver.session() as session:
                await upsert_pipeline_records(session, new_batch)

        # Forward ONLY new records (if not already present) to output queues
        for record in new_batch:
            unique_field = record.data.get_unique_model_field()
            unique_value = str(getattr(record.data, unique_field))
            context.added_nodes[type(record.data).__name__].add(unique_value)
            for output_queue in output_queues.values():
                await output_queue.put(record)

        # Update progress for all records
        if progress is not None and task_id is not None:
            progress.advance(task_id, len(batch))


class MilvusUpsertStage:
    """DEPRECATED: Batches PipelineRecord embeddings and upserts to Milvus.

    This stage uses the legacy VectorDB class. For new code, use:
    - pyeed.db.milvus.bulk_insert() for column-oriented inserts
    - pyeed.embed.embed_from_db.run_embedding_job() for full pipeline
    """

    def __init__(
        self,
        vector_db: Any,  # VectorDB (legacy)
        collection_name: str,
        batch_size: int = 2000,
    ):
        """Initialize Milvus sink.

        Args:
            vector_db: VectorDB instance (legacy)
            collection_name: Milvus collection name
            batch_size: Records to accumulate before writing
        """
        warnings.warn(
            "MilvusUpsertStage is deprecated. Use pyeed.db.milvus.bulk_insert()",
            DeprecationWarning,
            stacklevel=2,
        )
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
        batch: list[PipelineRecord[BaseNode]] = []

        while True:
            item = await input_queue.get()

            # Log the current size of the input queue
            logger.debug(f"MilvusUpsertStage input queue size: {input_queue.qsize()}")

            if progress is not None and task_id is not None:
                total = context.stats.get("total")
                if total is not None:
                    progress.update(task_id, total=total)

            if item is SENTINEL:
                # Flush remaining batch
                if batch:
                    await self._process_batch(batch, output_queues, progress, task_id)
                break

            # Skip records without embeddings
            if not item.embeddings:
                if progress is not None and task_id is not None:
                    progress.advance(task_id, 1)
                continue

            batch.append(item)

            if len(batch) >= self.batch_size:
                await self._process_batch(batch, output_queues, progress, task_id)
                batch = []

    async def _process_batch(
        self,
        batch: list[PipelineRecord[BaseNode]],
        output_queues: dict[str, asyncio.Queue[Any]],
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Process a batch: check existing protein_ids, filter, insert new ones, forward all."""
        if not batch:
            return

        # Extract protein_ids from batch
        protein_ids = [record.data.id for record in batch]

        # Check which protein_ids already exist in Milvus
        existing_ids: set[str] = set()
        with contextlib.suppress(Exception):
            # If check fails, proceed with all records (fail-safe)
            existing_ids = await asyncio.to_thread(
                self.vector_db._check_existing_ids,
                protein_ids,
            )

        # Filter batch: keep only records where protein_id doesn't exist
        new_batch: list[PipelineRecord[BaseNode]] = [
            record for record in batch if record.data.id not in existing_ids
        ]

        # Insert only new records
        if new_batch:
            embedding_batch = self._convert_to_embedding_batch(new_batch)
            await self.vector_db.insert_async(
                self.collection_name,
                embedding_batch,
                include_sequence=False,  # Don't store sequences
            )

        # Forward ALL records (both new and existing) to output queues
        for record in batch:
            for output_queue in output_queues.values():
                await output_queue.put(record)

        # Update progress for all records
        if progress is not None and task_id is not None:
            progress.advance(task_id, len(batch))

    def _convert_to_embedding_batch(
        self,
        records: list[PipelineRecord[BaseNode]],
    ) -> dict[str, Any]:
        """Convert PipelineRecords to embedding batch dict.

        Args:
            records: List of PipelineRecord with embeddings

        Returns:
            Dict with protein_ids, sequences, and embeddings by pooling method
        """
        batch: dict[str, Any] = {
            "protein_ids": [],
            "sequences": [],
            "embeddings": {},
        }

        for record in records:
            batch["protein_ids"].append(record.data.id)
            batch["sequences"].append(record.data.sequence)

            for pooling_name, embedding in record.embeddings.items():
                if pooling_name not in batch["embeddings"]:
                    batch["embeddings"][pooling_name] = []
                batch["embeddings"][pooling_name].append(embedding)

        return batch
