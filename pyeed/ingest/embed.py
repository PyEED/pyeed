"""Embedding computation and storage functions for proteins."""

from __future__ import annotations

import asyncio
import re

from loguru import logger
from rich.progress import Progress, TaskID

from ..db.neo4j import Database
from ..embedding.esm2_async import ESM2Embedder
from ..embedding.pooling import mean_pooling
from ..embedding.sinks import EmbeddingSink
from .progress import create_progress

__all__ = [
    "embed_proteins",
]


async def embed_proteins(
    db: Database,
    *,
    model_name: str = "facebook/esm2_t33_650M_UR50D",
    pooling_methods: list | None = None,
    sinks: list[EmbeddingSink] | None = None,
    batch_size: int = 16,
    chunk_size: int = 1000,
    model_dtype: str = "float16",
    return_dtype: str = "float16",
    show_progress: bool = True,
    progress: Progress | None = None,
    only_missing: bool = True,
) -> int:
    """Embed proteins and write to configured sinks.

    This function:
    1. Queries database for proteins needing embeddings
    2. Computes embeddings using ESM2Embedder with specified pooling methods
    3. Writes embeddings to all configured sinks (default: Neo4j as vector properties)

    Args:
        db: Database instance
        model_name: ESM2 model identifier (e.g., "facebook/esm2_t33_650M_UR50D")
        pooling_methods: List of pooling functions (default: [mean_pooling])
        sinks: List of EmbeddingSink instances (default)
        batch_size: Sequences per GPU batch
        chunk_size: Sequences per embedding chunk (controls memory usage)
        model_dtype: Model computation dtype (float32/float16/bfloat16)
        return_dtype: Returned embedding dtype (float32/float16/bfloat16)
        show_progress: Display progress bar
        progress: Existing Progress instance to share progress context
        only_missing: Only embed proteins without existing embeddings

    Returns:
        Number of proteins embedded
    """
    # Set defaults
    pooling_methods = pooling_methods or [mean_pooling]

    # Build property names for checking if embeddings exist
    # Following the pattern: vec__{model_name}__{pooling_method}

    model_slug = re.sub(r"[^a-z0-9_]", "_", model_name.lower())
    pool_names = [
        re.sub(r"[^a-z0-9_]", "_", getattr(p, "__name__", "custom").lower())
        for p in pooling_methods
    ]
    vec_props = [f"vec__{model_slug}__{pn}" for pn in pool_names]

    # Query for proteins to embed
    if only_missing:
        # Check if ANY of the vector properties are missing
        where_clauses = " OR ".join([f"p.{vp} IS NULL" for vp in vec_props])
        query = f"""
        MATCH (p:Protein)
        WHERE p.sequence IS NOT NULL
          AND ({where_clauses})
        RETURN p.sequence_id AS accession, p.sequence AS sequence
        """
        count_query = f"""
        MATCH (p:Protein)
        WHERE p.sequence IS NOT NULL
          AND ({where_clauses})
        RETURN count(p) AS n
        """
    else:
        query = """
        MATCH (p:Protein)
        WHERE p.sequence IS NOT NULL
        RETURN p.sequence_id AS accession, p.sequence AS sequence
        """
        count_query = """
        MATCH (p:Protein)
        WHERE p.sequence IS NOT NULL
        RETURN count(p) AS n
        """

    # Count total proteins
    count_result = db.query(count_query)
    total_proteins = count_result[0]["n"] if count_result else 0

    if total_proteins == 0:
        logger.info("No proteins to embed")
        return 0

    logger.info(
        f"Embedding {total_proteins} proteins with {len(pooling_methods)} pooling method(s)"
    )

    # Initialize embedder
    embedder = ESM2Embedder(
        model_name=model_name,
        pooling_methods=pooling_methods,
        model_dtype=model_dtype,
        return_dtype=return_dtype,
        normalize=True,
        verbose=False,  # We'll show our own progress bar
    )

    await embedder.initialize()

    # Connect all sinks
    await asyncio.gather(*[sink.connect() for sink in sinks])

    try:
        # Producer-consumer pattern with queue
        queue: asyncio.Queue[tuple[list[str], list[str]] | None] = asyncio.Queue(maxsize=chunk_size)

        async def producer(fetch_task_id: TaskID | None, prog: Progress | None) -> None:
            """Fetch proteins from database in chunks."""
            batch_seqs: list[str] = []
            batch_accs: list[str] = []

            async for rec in db.async_query_iter(query):
                batch_seqs.append(rec["sequence"])
                batch_accs.append(rec["accession"])

                if len(batch_seqs) >= chunk_size:
                    await queue.put((batch_seqs, batch_accs))
                    if fetch_task_id is not None and prog is not None:
                        prog.update(fetch_task_id, advance=len(batch_seqs))
                    batch_seqs = []
                    batch_accs = []

            # Flush remaining
            if batch_seqs:
                await queue.put((batch_seqs, batch_accs))
                if fetch_task_id is not None and prog is not None:
                    prog.update(fetch_task_id, advance=len(batch_seqs))

            await queue.put(None)  # Sentinel

        async def consumer(save_task_id: TaskID | None, prog: Progress | None) -> int:
            """Embed sequences and write to all sinks."""
            total_embedded = 0

            while True:
                item = await queue.get()
                if item is None:
                    break

                sequences, accessions = item

                # Compute embeddings for this chunk
                result = await embedder.embed_batch(sequences, accessions, batch_size=batch_size)

                # Prepare metadata
                metadata = {
                    "model_name": model_name,
                    "layer_index": -1,  # ESM2 uses last layer
                }

                # Write to all sinks in parallel
                await asyncio.gather(*[sink.write_batch(result, metadata) for sink in sinks])

                total_embedded += len(accessions)

                if save_task_id is not None and prog is not None:
                    prog.update(save_task_id, advance=len(accessions))

            return total_embedded

        # Execute pipeline with progress tracking
        if show_progress:
            progress_instance = create_progress(progress=progress)
            # Only use 'with' context if we created a new progress instance
            if progress is None:
                with progress_instance:
                    fetch_task = progress_instance.add_task("Fetch Proteins", total=total_proteins)
                    save_task = progress_instance.add_task("Embed & Save", total=total_proteins)
                    _, embedded = await asyncio.gather(
                        producer(fetch_task, progress_instance),
                        consumer(save_task, progress_instance),
                    )
            else:
                # Shared progress - don't use context manager
                fetch_task = progress_instance.add_task("Fetch Proteins", total=total_proteins)
                save_task = progress_instance.add_task("Embed & Save", total=total_proteins)
                _, embedded = await asyncio.gather(
                    producer(fetch_task, progress_instance),
                    consumer(save_task, progress_instance),
                )
        else:
            # No progress tracking
            _, embedded = await asyncio.gather(
                producer(None, None),
                consumer(None, None),
            )

        logger.info(f"Embedded {embedded} proteins")
        return embedded

    finally:
        # Cleanup
        await asyncio.gather(*[sink.disconnect() for sink in sinks])
        await embedder.cleanup()
