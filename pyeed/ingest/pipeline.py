from __future__ import annotations

import asyncio
from collections.abc import AsyncIterator, Callable

import numpy as np
from loguru import logger
from pymilvus.exceptions import DescribeCollectionException

from ..db.milvus import VectorDB
from ..embed.esm2 import ESM2Embedder
from ..embed.types import NP_DTYPE_MAP, EmbeddingBatch
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
        max_reader_queue: int = 256,
        max_writer_queue: int = 256,
        include_sequence: bool = True,
        header_extractor: Callable[[str], str] | None = None,
    ):
        """Initialize embedding pipeline.

        Args:
            embedder: ESM2Embedder instance (must be initialized)
            vector_db: VectorDB instance for Milvus operations
            collection_name: Name of Milvus collection to write to
            chunk_size: FASTA entries per chunk (default: 3200)
            batch_size: Sequences per GPU batch (default: 32)
            write_batch_size: Batches to accumulate before Milvus write (default: 10)
            max_reader_queue: Max batches buffered after reading (default: 256)
            max_writer_queue: Max batches buffered after embedding (default: 256)
            include_sequence: Whether to include sequence in Milvus records
            header_extractor: Optional function to extract protein_id from FASTA header
                             Signature: (header: str) -> str
                             If None, uses entire header as protein_id
        """
        self.embedder = embedder
        self.vector_db = vector_db
        self.collection_name = collection_name
        self.chunk_size = chunk_size
        self.batch_size = batch_size
        self.write_batch_size = write_batch_size
        self.include_sequence = include_sequence
        self.header_extractor = header_extractor

        # Two queues only: reader → embedder → writer
        self.reader_queue: asyncio.Queue[tuple[list[str], list[str]] | object] = asyncio.Queue(
            maxsize=max_reader_queue
        )
        self.writer_queue: asyncio.Queue[EmbeddingBatch | object] = asyncio.Queue(
            maxsize=max_writer_queue
        )

    async def _reader_worker(
        self,
        fasta_path: str,
        read_progress_reporter: ProgressReporter,
        embed_progress_reporter: ProgressReporter,
        write_progress_reporter: ProgressReporter,
    ) -> None:
        """Read FASTA, check existing IDs, sort by length, create batches, emit to queue.

        Splits sequences by embedder's max_length and filters existing IDs:
        - Existing IDs → skipped (not processed)
        - OK sequences (≤ max_length) → reader_queue → embedder
        - Long sequences (> max_length) → writer_queue with zero embeddings

        Args:
            fasta_path: Path to FASTA file
            read_progress_reporter: Progress reporter for reading
            embed_progress_reporter: Progress reporter for embedding
            write_progress_reporter: Progress reporter for writing
        """
        max_length = self.embedder.max_length

        async for chunk_dict in read_fasta_chunks_async(
            fasta_path,
            chunk_size=self.chunk_size,
            header_extractor=self.header_extractor,
        ):
            read_progress_reporter(advance=len(chunk_dict))

            # Collect all protein IDs from chunk
            all_protein_ids = list(chunk_dict.keys())

            # Check for existing IDs in Milvus (run in thread pool)
            try:
                existing_ids = await asyncio.to_thread(
                    self.vector_db._check_existing_ids,
                    all_protein_ids,
                )
            except DescribeCollectionException as e:
                logger.error(f"Error checking existing IDs: {e}")
                existing_ids = set()

            # Filter out existing IDs
            new_pairs: list[tuple[str, str]] = [
                (pid, seq) for pid, seq in chunk_dict.items() if pid not in existing_ids
            ]

            # Count skipped sequences for progress adjustment
            skipped_count = len(existing_ids)
            if skipped_count > 0:
                logger.debug(f"Skipping {skipped_count} existing sequences in chunk")
                # Adjust embed progress total
                total_sequences = embed_progress_reporter.progress._tasks[
                    embed_progress_reporter.task_id
                ].total
                new_total = total_sequences - skipped_count
                embed_progress_reporter.progress.update(
                    embed_progress_reporter.task_id, total=new_total
                )
                # Adjust write progress total
                write_total = write_progress_reporter.progress._tasks[
                    write_progress_reporter.task_id
                ].total
                write_new_total = write_total - skipped_count
                write_progress_reporter.progress.update(
                    write_progress_reporter.task_id, total=write_new_total
                )

            # Split new sequences by length
            ok_pairs: list[tuple[str, str]] = []
            long_pairs: list[tuple[str, str]] = []

            for protein_id, sequence in new_pairs:
                if len(sequence) <= max_length:
                    ok_pairs.append((protein_id, sequence))
                else:
                    long_pairs.append((protein_id, sequence))

            # Process OK sequences (normal flow)
            if ok_pairs:
                # Sort by length in thread pool (descending - longest first)
                sorted_pairs = await asyncio.to_thread(
                    sorted,
                    ok_pairs,
                    key=lambda x: len(x[1]),
                    reverse=True,
                )

                # Create batches
                for i in range(0, len(sorted_pairs), self.batch_size):
                    batch = sorted_pairs[i : i + self.batch_size]
                    accessions = [acc for acc, _ in batch]
                    sequences = [seq for _, seq in batch]
                    await self.reader_queue.put((sequences, accessions))

            # Process long sequences (skip embedder, use zero embeddings)
            if long_pairs:
                # Create batches of long sequences
                for i in range(0, len(long_pairs), self.batch_size):
                    batch = long_pairs[i : i + self.batch_size]
                    protein_ids = [pid for pid, _ in batch]
                    sequences = [seq for _, seq in batch]

                    # Create EmbeddingBatch with zero embeddings
                    zero_batch = self._create_zero_embedding_batch(protein_ids, sequences)
                    await self.writer_queue.put(zero_batch)
                    # Adjust embed progress for long sequences (skip embedding)
                    total_sequences = embed_progress_reporter.progress._tasks[
                        embed_progress_reporter.task_id
                    ].total
                    new_total = total_sequences - len(long_pairs)
                    embed_progress_reporter.progress.update(
                        embed_progress_reporter.task_id, total=new_total
                    )

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

    def _create_zero_embedding_batch(
        self,
        protein_ids: list[str],
        sequences: list[str],
    ) -> EmbeddingBatch:
        """Create EmbeddingBatch with zero embeddings for long sequences.

        Args:
            protein_ids: List of protein IDs
            sequences: List of sequences

        Returns:
            EmbeddingBatch with zero embeddings (single vector per sequence)
        """
        batch = EmbeddingBatch()
        batch.protein_ids = protein_ids
        batch.sequences = sequences

        # Get embedding dimension from embedder (model hidden_size)
        embedding_dim = self.embedder.embedding_dim

        # Get return dtype
        dtype = NP_DTYPE_MAP[self.embedder.return_dtype]

        # Create zero vector for each sequence (just model output dimension)
        # Use first pooling method name or "zero" as key
        pool_name = self.embedder.pooling_configs[0][0] if self.embedder.pooling_configs else "zero"

        # Create zero embeddings: list of zero vectors, one per sequence
        batch.embeddings[pool_name] = [np.zeros(embedding_dim, dtype=dtype) for _ in protein_ids]

        print("Created zero embedding batch")
        return batch

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
                self._reader_worker(
                    fasta_path,
                    ProgressReporter(progress, read_task),
                    ProgressReporter(progress, embed_task),
                    ProgressReporter(progress, write_task),
                ),
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

    def extract_uniprot_id(s: str) -> str:
        """Extract UniProt ID from header using regex."""
        return s.split("|")[1]

    async def main() -> None:
        embedder = ESM2Embedder(
            model_name="facebook/esm2_t33_650M_UR50D",
            model_dtype="float32",
            return_dtype="float32",
            n_gpus=2,
            pooling_methods=[mean_pooling],
        )
        await embedder.initialize()

        vector_db = VectorDB(uri="http://localhost:19530", token="root:Milvus")
        pipeline = EmbeddingPipeline(
            embedder,
            vector_db,
            "test",
            batch_size=16,
            chunk_size=16 * 400,
            write_batch_size=16,
            header_extractor=extract_uniprot_id,
        )

        await pipeline.run(
            "/home/mha/projects/proteingraph/downloads/uniprot_sprot.fasta",
        )

        await embedder.cleanup()

    asyncio.run(main())
