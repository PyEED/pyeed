"""Three-stage async embedding pipeline for FASTA → Embedding → Milvus."""

from __future__ import annotations

import asyncio
import re
from collections.abc import AsyncIterator

from loguru import logger

from ..db.milvus import VectorDB
from ..embedding.esm2_async import ESM2Embedder
from ..embedding.types import EmbeddingBatch
from ..utils import build_header_index, read_fasta_chunks_async
from .progress import ProgressReporter, create_progress

__all__ = [
    "EmbeddingPipeline",
]


class EmbeddingPipeline:
    """Three-stage async pipeline: FASTA Reading → Embedding → Milvus Writing.

    This design ensures maximum throughput by overlapping I/O and compute:
    - Reader reads FASTA, sorts by length, and creates batches
    - Embedder processes batches on GPU(s) as they arrive
    - Writer accumulates batches and writes to Milvus in groups

    Key insight: Use bounded queues to prevent memory explosion while
    maintaining concurrency.
    """

    SENTINEL = object()

    def __init__(  # noqa: PLR0913
        self,
        embedder: ESM2Embedder,
        vector_db: VectorDB,
        collection_name: str,
        *,
        chunk_size: int = 32 * 100,
        batch_size: int = 32,
        write_batch_size: int = 10,
        max_reader_queue: int = 2,
        max_writer_queue: int = 2,
        include_sequence: bool = True,
        header_pattern: str | re.Pattern[str] | None = None,
    ):
        """Initialize embedding pipeline.

        Args:
            embedder: ESM2Embedder instance (must be initialized)
            vector_db: VectorDB instance for Milvus operations
            collection_name: Name of Milvus collection to write to
            chunk_size: FASTA entries per chunk (default: 3200)
            batch_size: Sequences per GPU batch (default: 32)
            write_batch_size: Batches to accumulate before Milvus write (default: 10)
            max_reader_queue: Max batches buffered after reading (default: 2)
            max_writer_queue: Max batches buffered after embedding (default: 2)
            include_sequence: Whether to include sequence in Milvus records
            header_pattern: Optional regex pattern to extract protein_id from FASTA header
        """
        self.embedder = embedder
        self.vector_db = vector_db
        self.collection_name = collection_name
        self.chunk_size = chunk_size
        self.batch_size = batch_size
        self.write_batch_size = write_batch_size
        self.include_sequence = include_sequence
        self.header_pattern = header_pattern

        # Two queues only: reader → embedder → writer
        self.reader_queue: asyncio.Queue[tuple[list[str], list[str]] | object] = asyncio.Queue(
            maxsize=max_reader_queue
        )
        self.writer_queue: asyncio.Queue[EmbeddingBatch | object] = asyncio.Queue(
            maxsize=max_writer_queue
        )

        # Statistics
        self._total_written = 0

    async def _reader_worker(
        self,
        fasta_path: str,
        progress_reporter: ProgressReporter,
    ) -> None:
        """Read FASTA, sort by length, create batches, emit to queue.

        Args:
            fasta_path: Path to FASTA file
            progress_reporter: Progress reporter for reading
        """
        async for chunk_dict in read_fasta_chunks_async(
            fasta_path,
            chunk_size=self.chunk_size,
            header_pattern=self.header_pattern,
        ):
            # Sort by length in thread pool (descending - longest first)
            sorted_pairs = await asyncio.to_thread(
                sorted,
                list(chunk_dict.items()),
                key=lambda x: len(x[1]),
                reverse=True,
            )

            # Create batches
            for i in range(0, len(sorted_pairs), self.batch_size):
                batch = sorted_pairs[i : i + self.batch_size]
                accessions = [acc for acc, _ in batch]
                sequences = [seq for _, seq in batch]
                await self.reader_queue.put((sequences, accessions))
                progress_reporter(advance=len(sequences))

        await self.reader_queue.put(self.SENTINEL)
        logger.debug("Reader worker finished")

    async def _embedder_worker(
        self,
        progress_reporter: ProgressReporter,
    ) -> None:
        """Convert queue to async iterator, call embedder.embed_stream().

        Args:
            progress_reporter: Progress reporter for embedding
        """

        async def batch_iterator() -> AsyncIterator[tuple[list[str], list[str]]]:
            """Convert queue to async iterator."""
            while True:
                item = await self.reader_queue.get()
                if item is self.SENTINEL:
                    break
                yield item  # type: ignore

        async for embedding_batch in self.embedder.embed_stream(batch_iterator()):
            await self.writer_queue.put(embedding_batch)
            progress_reporter(advance=len(embedding_batch.protein_ids))

        await self.writer_queue.put(self.SENTINEL)
        logger.debug("Embedder worker finished")

    async def _writer_worker(
        self,
        progress_reporter: ProgressReporter,
    ) -> None:
        """Accumulate batches, write to Milvus in groups.

        Args:
            progress_reporter: Progress reporter for writing
        """
        buffer: list[EmbeddingBatch] = []

        while True:
            item = await self.writer_queue.get()

            if item is self.SENTINEL:
                # Flush remaining batches
                if buffer:
                    combined = self._combine_batches(buffer)
                    inserted = await self.vector_db.insert_async(
                        self.collection_name,
                        combined,
                        include_sequence=self.include_sequence,
                    )
                    self._total_written += inserted
                    progress_reporter(advance=inserted)
                break

            buffer.append(item)  # type: ignore

            # Flush when buffer is full
            if len(buffer) >= self.write_batch_size:
                combined = self._combine_batches(buffer)
                inserted = await self.vector_db.insert_async(
                    self.collection_name,
                    combined,
                    include_sequence=self.include_sequence,
                )
                self._total_written += inserted
                progress_reporter(advance=inserted)
                buffer.clear()

        logger.debug("Writer worker finished")

    def _combine_batches(self, batches: list[EmbeddingBatch]) -> EmbeddingBatch:
        """Combine multiple EmbeddingBatch into one.

        Args:
            batches: List of EmbeddingBatch to combine

        Returns:
            Single combined EmbeddingBatch
        """
        combined = EmbeddingBatch()
        for batch in batches:
            combined.protein_ids.extend(batch.protein_ids)
            combined.sequences.extend(batch.sequences)
            for pool_name, embeddings in batch.embeddings.items():
                if pool_name not in combined.embeddings:
                    combined.embeddings[pool_name] = []
                combined.embeddings[pool_name].extend(embeddings)
        return combined

    async def run(self, fasta_path: str) -> int:
        """Run the three-stage pipeline.

        Args:
            fasta_path: Path to FASTA file

        Returns:
            Total number of proteins embedded and written

        Raises:
            FileNotFoundError: If FASTA file doesn't exist
            RuntimeError: If pipeline errors occur
        """
        logger.info(
            f"Starting embedding pipeline for '{fasta_path}'",
            extra={
                "collection_name": self.collection_name,
                "chunk_size": self.chunk_size,
                "batch_size": self.batch_size,
                "write_batch_size": self.write_batch_size,
            },
        )

        # Reset statistics
        self._total_written = 0

        # Build index for progress tracking
        logger.info("Building header index...")
        offsets = await asyncio.to_thread(build_header_index, fasta_path)
        total_sequences = len(offsets)
        logger.info(f"Found {total_sequences} sequences in FASTA file")

        # Set up progress tracking
        with create_progress() as progress:
            read_task = progress.add_task("Read & Batch", total=total_sequences)
            embed_task = progress.add_task("Embed", total=total_sequences)
            write_task = progress.add_task("Write", total=total_sequences)

            await asyncio.gather(
                self._reader_worker(fasta_path, ProgressReporter(progress, read_task)),
                self._embedder_worker(ProgressReporter(progress, embed_task)),
                self._writer_worker(ProgressReporter(progress, write_task)),
            )

        logger.info(
            "Pipeline completed",
            extra={"total_written": self._total_written},
        )

        return self._total_written


if __name__ == "__main__":
    from ..db.milvus import VectorDB
    from ..embedding.pooling import mean_pooling

    async def main() -> None:
        embedder = ESM2Embedder(
            model_name="facebook/esm2_t33_650M_UR50D",
            model_dtype="float32",
            return_dtype="float32",
            pooling_methods=[mean_pooling],
        )
        await embedder.initialize()

        vector_db = VectorDB(uri="http://localhost:19530", token="root:Milvus")
        pipeline = EmbeddingPipeline(
            embedder,
            vector_db,
            "test",
            batch_size=24,
            chunk_size=24 * 100,
            write_batch_size=10,
        )

        await pipeline.run(
            "/home/mha/projects/proteingraph/downloads/uniprot_sprot.fasta",
        )

        await embedder.cleanup()

    asyncio.run(main())
