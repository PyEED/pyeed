from __future__ import annotations

import asyncio
from collections.abc import Callable
from typing import Any

from rich.progress import Progress, TaskID

from pyeed.ingest.core.pipeline import PipelineRecord

from ..core.protocol import SENTINEL, PipelineContext
from ..model import Protein
from ..sources.fasta import read_fasta_chunks_async


class FASTAReaderStage:
    """Reads FASTA file and emits Protein objects to queue."""

    def __init__(
        self,
        fasta_path: str,
        offsets: list[int],
        chunk_size: int,
        header_extractor: Callable[[str], str] | None,
    ):
        """Initialize FASTA reader.

        Args:
            fasta_path: Path to FASTA file
            chunk_size: Number of sequences per chunk
            header_extractor: Optional function to extract protein_id from header
        """
        self.fasta_path = fasta_path
        self.chunk_size = chunk_size
        self.header_extractor = header_extractor
        self.offsets = offsets

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Read FASTA and emit protein records."""

        async for chunk in read_fasta_chunks_async(
            self.fasta_path,
            offsets=self.offsets,
            chunk_size=self.chunk_size,
            header_extractor=self.header_extractor,
        ):
            if progress is not None and task_id is not None:
                progress.advance(task_id, len(chunk))

            for seq_id, sequence in chunk.items():
                protein = Protein(
                    id=seq_id,
                    sequence=sequence,
                    seq_length=len(sequence),
                )
                record = PipelineRecord(data=protein)
                for output_queue in output_queues.values():
                    await output_queue.put(record)

        for output_queue in output_queues.values():
            await output_queue.put(SENTINEL)
