"""Neo4j sink stage for batch processing PipelineRecords."""

from __future__ import annotations

import asyncio
from collections import defaultdict
from typing import Any

from rich.progress import Progress, TaskID

from pyeed.db.neo4j import GraphDB
from pyeed.db.queries import upsert_pipeline_records
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
        """Consume PipelineRecords, batch, and upsert to Neo4j with relationships."""
        input_queue: asyncio.Queue[PipelineRecord[PyeedBase]] = next(iter(input_queues.values()))

        item_dict: dict[str, set[str]] = defaultdict(set)
        batch: list[PipelineRecord[PyeedBase]] = []

        while True:
            item = await input_queue.get()

            if item is SENTINEL:
                # Flush remaining batch
                if batch:
                    async with self.db.async_driver.session() as session:
                        await upsert_pipeline_records(session, batch)
                    # Track all parent nodes
                    for record in batch:
                        unique_field = record.data.get_unique_model_field()
                        unique_value = str(getattr(record.data, unique_field))
                        context.added_nodes[type(record.data).__name__].add(unique_value)
                    if progress and task_id:
                        progress.advance(task_id, len(batch))
                break

            unique_field = item.data.get_unique_model_field()
            unique_value = str(getattr(item.data, unique_field))
            item_dict[type(item.data).__name__].add(unique_value)

            batch.append(item)

            if len(batch) >= self.batch_size:
                async with self.db.async_driver.session() as session:
                    await upsert_pipeline_records(session, batch)
                # Track all parent nodes
                for record in batch:
                    unique_field = record.data.get_unique_model_field()
                    unique_value = str(getattr(record.data, unique_field))
                    context.added_nodes[type(record.data).__name__].add(unique_value)
                if progress and task_id:
                    progress.advance(task_id, len(batch))
                batch = []
