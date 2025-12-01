from __future__ import annotations

import asyncio
import re
from collections.abc import Callable
from typing import Any

from loguru import logger
from rich.progress import Progress, TaskID

from pyeed.ingest.core.pipeline import IngestItem

from ..core.protocol import SENTINEL, PipelineContext
from ..model import Protein
from ..sources.fasta import read_fasta_chunks_async


class FASTAReaderStage:
    """Reads FASTA file and emits Protein objects to queue."""

    VALID_PROTEIN_ID_PATTERN = re.compile(r"^[a-zA-Z0-9_|.\-]+$")
    MAX_PROTEIN_ID_LENGTH = 255

    def __init__(
        self,
        fasta_path: str,
        offsets: list[int],
        chunk_size: int,
        header_extractor: Callable[[str], str] | None,
        taxon_extractor: Callable[[str], str] | None = None,
    ) -> None:
        """Initialize FASTA reader.

        Args:
            fasta_path: Path to FASTA file
            chunk_size: Number of sequences per chunk
            header_extractor: Optional function to extract protein_id from header
            taxon_extractor: Optional function to extract taxon_id from header
        """
        self.fasta_path = fasta_path
        self.chunk_size = chunk_size
        self.header_extractor = header_extractor
        self.taxon_extractor = taxon_extractor
        self.offsets = offsets

    def _is_valid_protein_id(self, protein_id: str) -> bool:
        """Validate protein ID format.

        Args:
            protein_id: Protein ID to validate

        Returns:
            True if valid, False otherwise
        """
        return (
            bool(self.VALID_PROTEIN_ID_PATTERN.match(protein_id))
            and len(protein_id) <= self.MAX_PROTEIN_ID_LENGTH
        )

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Read FASTA and emit protein records.

        Args:
            input_queues: Input queues (unused for reader stage)
            output_queues: Output queues to emit PipelineRecords
            context: Pipeline context for tracking state
            progress: Rich progress instance for updates
            task_id: Task ID for progress tracking
        """

        skipped_count = 0

        async for chunk in read_fasta_chunks_async(
            self.fasta_path,
            offsets=self.offsets,
            chunk_size=self.chunk_size,
            header_extractor=self.header_extractor,
            taxon_extractor=self.taxon_extractor,
        ):
            valid_count = 0
            for seq_id, sequence, taxon_id in chunk:
                # Validate protein ID before creating Protein object
                if not self._is_valid_protein_id(seq_id):
                    logger.error(
                        f"Invalid protein ID format, skipping: {seq_id[:100]!r} "
                        f"(contains illegal characters or exceeds max length)"
                    )
                    skipped_count += 1
                    continue

                protein = Protein(
                    id=seq_id,
                    sequence=sequence,
                    seq_length=len(sequence),
                )  # type: ignore

                relations = {"TAXON": [taxon_id]} if taxon_id else None
                record = IngestItem(node=protein, relations=relations)
                for output_queue in output_queues.values():
                    await output_queue.put(record)

                valid_count += 1

            if progress is not None and task_id is not None:
                progress.advance(task_id, valid_count)

        if skipped_count > 0:
            logger.warning(
                f"Skipped {skipped_count} records with invalid protein IDs during FASTA reading"
            )

        for output_queue in output_queues.values():
            await output_queue.put(SENTINEL)
