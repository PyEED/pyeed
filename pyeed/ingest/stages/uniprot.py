"""UniProt and InterPro reader stages for fetching protein data from REST APIs."""

from __future__ import annotations

import asyncio
from typing import Any

import httpx
from loguru import logger
from rich.progress import Progress, TaskID

from pyeed.ingest.core.protocol import SENTINEL, PipelineContext
from pyeed.ingest.sources.uniprot import UniProtAdapter


class UniProtReaderStage:
    """Fetches proteins from UniProt REST API and emits PipelineRecords.

    Fetches protein data by accession IDs, including rich metadata such as
    GO annotations, sequence annotations, EC numbers, and Rhea reaction IDs.
    """

    def __init__(
        self,
        accessions: list[str],
        chunk_size: int = 50,
        size_per_page: int = 50,
    ) -> None:
        """Initialize UniProt reader.

        Args:
            accessions: List of UniProt accession IDs to fetch
            chunk_size: Number of accessions per API query batch (default: 50)
            size_per_page: Results per page for pagination (default: 50)
        """
        self.accessions = accessions
        self.chunk_size = chunk_size
        self.size_per_page = size_per_page
        self.adapter = UniProtAdapter()

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Fetch proteins from UniProt and emit PipelineRecords.

        Args:
            input_queues: Input queues (unused for reader stage)
            output_queues: Output queues to emit PipelineRecords
            context: Pipeline context for tracking state
            progress: Rich progress instance for updates
            task_id: Task ID for progress tracking
        """
        processed_count = 0
        failed_count = 0

        async with httpx.AsyncClient() as client:
            try:
                async for rec in self.adapter.fetch_accessions(
                    client=client,
                    accessions=self.accessions,
                    chunk_size=self.chunk_size,
                    size_per_page=self.size_per_page,
                ):
                    try:
                        # Map UniProt record to PipelineRecord with children
                        pipeline_record = self.adapter.map(rec)

                        # Emit to all output queues
                        for output_queue in output_queues.values():
                            await output_queue.put(pipeline_record)

                        processed_count += 1

                        if progress is not None and task_id is not None:
                            progress.advance(task_id, 1)

                    except Exception as e:
                        accession = rec.get("primaryAccession", "unknown")
                        logger.error(
                            f"Failed to process UniProt record {accession}: {e}",
                            exc_info=True,
                        )
                        failed_count += 1

            except Exception as e:
                logger.error(f"UniProt fetch failed: {e}", exc_info=True)

        if failed_count > 0:
            logger.warning(
                f"UniProt reader completed: {processed_count} processed, {failed_count} failed"
            )
        else:
            logger.info(f"UniProt reader completed: {processed_count} proteins fetched")

        # Send SENTINEL to signal completion
        for output_queue in output_queues.values():
            await output_queue.put(SENTINEL)


class InterProReaderStage:
    """Fetches all proteins for an InterPro ID and emits PipelineRecords.

    First queries the UniProt SPARQL endpoint to get all accession IDs
    associated with an InterPro family, then fetches full protein data
    for each accession.
    """

    def __init__(
        self,
        interpro_id: str,
        chunk_size: int = 50,
        size_per_page: int = 50,
        distinct: bool = True,
        limit: int = 1_000_000_000,
    ) -> None:
        """Initialize InterPro reader.

        Args:
            interpro_id: InterPro ID (e.g., "IPR002133")
            chunk_size: Number of accessions per API query batch (default: 50)
            size_per_page: Results per page for pagination (default: 50)
            distinct: Use SELECT DISTINCT in SPARQL query (default: True)
            limit: Maximum number of accessions to fetch (default: 1 billion)
        """
        self.interpro_id = interpro_id
        self.chunk_size = chunk_size
        self.size_per_page = size_per_page
        self.distinct = distinct
        self.limit = limit
        self.adapter = UniProtAdapter()

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Fetch proteins from InterPro and emit PipelineRecords.

        Args:
            input_queues: Input queues (unused for reader stage)
            output_queues: Output queues to emit PipelineRecords
            context: Pipeline context for tracking state
            progress: Rich progress instance for updates
            task_id: Task ID for progress tracking
        """
        async with httpx.AsyncClient() as client:
            # Step 1: Get all accessions for this InterPro ID
            logger.info(f"Fetching accessions for {self.interpro_id}")
            accessions = await self.adapter.fetch_accessions_by_interpro(
                client=client,
                ipr=self.interpro_id,
                distinct=self.distinct,
                limit=self.limit,
            )
            logger.info(f"Found {len(accessions)} accessions for {self.interpro_id}")

            # Update progress total if available
            if progress is not None and task_id is not None:
                progress.update(task_id, total=len(accessions))

            # Step 2: Fetch protein details
            processed_count = 0
            failed_count = 0

            try:
                async for rec in self.adapter.fetch_accessions(
                    client=client,
                    accessions=accessions,
                    chunk_size=self.chunk_size,
                    size_per_page=self.size_per_page,
                ):
                    try:
                        # Map UniProt record to PipelineRecord with children
                        pipeline_record = self.adapter.map(rec)

                        # Emit to all output queues
                        for output_queue in output_queues.values():
                            await output_queue.put(pipeline_record)

                        processed_count += 1

                        if progress is not None and task_id is not None:
                            progress.advance(task_id, 1)

                    except Exception as e:
                        accession = rec.get("primaryAccession", "unknown")
                        logger.error(
                            f"Failed to process UniProt record {accession}: {e}",
                            exc_info=True,
                        )
                        failed_count += 1

            except Exception as e:
                logger.error(f"InterPro fetch failed: {e}", exc_info=True)

        if failed_count > 0:
            logger.warning(
                f"InterPro reader completed: {processed_count} processed, {failed_count} failed"
            )
        else:
            logger.info(
                f"InterPro reader completed: {processed_count} proteins fetched "
                f"for {self.interpro_id}"
            )

        # Send SENTINEL to signal completion
        for output_queue in output_queues.values():
            await output_queue.put(SENTINEL)
