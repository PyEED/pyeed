import asyncio
import gc
import logging

import numpy as np
import torch
from numpy.typing import NDArray
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TaskID,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from pyeed.database import Database
from pyeed.embeddings import ESM2Embedder
from pyeed.model import Protein
from pyeed.model.embedding import Embedding

logger = logging.getLogger(__name__)


def mean_pooling(
    hidden_states: torch.Tensor, attention_mask: torch.Tensor | None = None
) -> torch.Tensor:
    """
    Apply mean pooling to hidden states.

    Args:
        hidden_states: Shape [batch, seq_len, hidden_dim]
        attention_mask: Shape [batch, seq_len], 1 for real tokens, 0 for padding

    Returns:
        Pooled embeddings [batch, hidden_dim]
    """
    if attention_mask is None:
        # Simple mean over sequence dimension
        return hidden_states.mean(dim=1)

    # Expand attention mask to match hidden_states dimensions
    # attention_mask: [batch, seq_len] -> [batch, seq_len, 1]
    expanded_mask = attention_mask.unsqueeze(-1).expand(hidden_states.size()).float()

    # Sum embeddings, weighted by mask
    sum_embeddings = torch.sum(hidden_states * expanded_mask, dim=1)

    # Sum mask values to get count of real tokens
    sum_mask = torch.clamp(expanded_mask.sum(dim=1), min=1e-9)

    # Average
    return sum_embeddings / sum_mask


def normalize_embeddings(embeddings: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    L2-normalize embeddings.

    Args:
        embeddings: Shape [batch, dim] or [dim]

    Returns:
        Normalized embeddings with same shape
    """
    if embeddings.ndim == 1:
        embeddings = embeddings.reshape(1, -1)
        squeeze = True
    else:
        squeeze = False

    norm = np.linalg.norm(embeddings, axis=1, keepdims=True)
    norm[norm == 0] = 1.0  # Avoid division by zero
    normalized = embeddings / norm

    if squeeze:
        return np.asarray(normalized[0], dtype=np.float64)
    return np.asarray(normalized, dtype=np.float64)


async def get_proteins_without_embeddings(
    db: Database,
    model_name: str,
    pooling_method: str = "mean_pooling",
    layer_index: int = -1,
) -> list[tuple[str, str]]:
    """
    Query proteins that don't have embeddings for specified model/pooling/layer.

    Returns:
        List of (sequence_id, sequence) tuples
    """
    # Build the vector property name to check
    model_name_clean = model_name.replace("/", "_").replace("-", "_").lower()
    pooling_clean = pooling_method.lower()
    vec_prop = f"vec__{model_name_clean}__{pooling_clean}"

    query = """
    MATCH (p:Protein)
    WHERE p.sequence IS NOT NULL
      AND NOT EXISTS {
        MATCH (p)-[:HAS_EMBEDDING]->(e:Embedding)
        WHERE e.model_name = $model_name
          AND e.pooling_method = $pooling_method
          AND e.layer_index = $layer_index
      }
    RETURN p.sequence_id AS sequence_id, p.sequence AS sequence
    """

    logger.info(
        "Querying proteins without embeddings",
        extra={
            "model_name": model_name,
            "pooling_method": pooling_method,
            "layer_index": layer_index,
        },
    )

    results = []
    async for record in db.async_query_iter(
        query,
        model_name=model_name_clean,
        pooling_method=pooling_clean,
        layer_index=layer_index,
    ):
        results.append((record["sequence_id"], record["sequence"]))

    logger.info("Found proteins without embeddings", extra={"count": len(results)})
    return results


async def embed_proteins(
    db: Database,
    model_name: str = "facebook/esm2_t33_650M_UR50D",
    pooling_method: str = "mean_pooling",
    layer_index: int = -1,
    gpu_batch_size: int = 8,
    stream_chunk_size: int = 1000,
    db_batch_size: int = 100,
    device_ids: list[int] | None = None,
) -> None:
    """
    Embed all proteins that don't have embeddings yet.

    Uses async streaming: as soon as a batch is computed, it's saved to DB.

    Args:
        db: Database instance
        model_name: ESM-2 model identifier
        pooling_method: Pooling function name
        layer_index: Layer to extract (-1 = last)
        gpu_batch_size: Sequences processed per GPU batch (memory-dependent)
        stream_chunk_size: Sequences per stream chunk (affects progress granularity)
        db_batch_size: Number of proteins accumulated before DB write
        device_ids: Specific GPU IDs to use (None = all available)
    """
    # Get proteins that need embeddings
    proteins_todo = await get_proteins_without_embeddings(
        db, model_name, pooling_method, layer_index
    )

    if not proteins_todo:
        logger.info("No proteins to embed")
        return

    accessions = [acc for acc, _ in proteins_todo]
    sequences = [seq for _, seq in proteins_todo]

    logger.info(
        "Starting embedding pipeline",
        extra={
            "num_proteins": len(proteins_todo),
            "gpu_batch_size": gpu_batch_size,
            "stream_chunk_size": stream_chunk_size,
            "db_batch_size": db_batch_size,
        },
    )

    # Limit queue size to prevent memory explosion - only hold ~2 DB batches worth
    # This creates backpressure so GPU doesn't run too far ahead of DB writes
    max_queue_items = min(db_batch_size * 2, 200)
    queue: asyncio.Queue[tuple[str, Embedding] | None] = asyncio.Queue(maxsize=max_queue_items)

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        transient=False,
    )

    async def producer(task_id: TaskID) -> None:
        """Generate embeddings and stream to queue."""
        async with ESM2Embedder(
            model_name=model_name,
            pooling=mean_pooling,
            layer_index=layer_index,
            device_ids=device_ids,
        ) as embedder:
            # Track position to map embeddings back to accessions
            chunk_start = 0
            async for embeddings_chunk in embedder.embed_stream(
                sequences, accessions, gpu_batch_size, stream_chunk_size
            ):
                # Put embeddings in queue (this blocks if queue is full - backpressure!)
                chunk_size = len(embeddings_chunk)
                for i, emb in enumerate(embeddings_chunk):
                    acc = accessions[chunk_start + i]
                    await queue.put((acc, emb))

                # Update progress once per chunk
                progress.update(task_id, advance=chunk_size)
                progress.refresh()
                chunk_start += chunk_size

                # Force garbage collection after each chunk to prevent memory buildup
                # This is critical when processing large datasets
                del embeddings_chunk
                gc.collect()
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

        # Signal done
        await queue.put(None)

    async def consumer(task_id: TaskID) -> None:
        """Attach embeddings to existing proteins in database."""
        # Build batches: sequence_id -> [Embedding, ...]
        batch: dict[str, list[Embedding]] = {}

        while True:
            item = await queue.get()
            if item is None:
                # Final flush
                if batch:
                    await db.attach(Protein, batch)
                    progress.update(task_id, advance=sum(len(v) for v in batch.values()))
                    batch.clear()
                    gc.collect()  # Final cleanup
                break

            # Unpack (accession, embedding) tuple
            sequence_id, embedding = item

            # Group embeddings by protein accession
            if sequence_id not in batch:
                batch[sequence_id] = []
            batch[sequence_id].append(embedding)

            # Save batch when it reaches target size
            if len(batch) >= db_batch_size:
                await db.attach(Protein, batch)
                count = sum(len(v) for v in batch.values())
                progress.update(task_id, advance=count)
                progress.refresh()
                batch.clear()
                # Aggressive garbage collection to free memory ASAP
                gc.collect()

    with progress:
        embed_task = progress.add_task("[cyan]Compute Embeddings", total=len(proteins_todo))
        save_task = progress.add_task("[green]Save to DB", total=len(proteins_todo))

        await asyncio.gather(
            producer(embed_task),
            consumer(save_task),
        )

    logger.info("Embedding pipeline complete", extra={"num_proteins": len(proteins_todo)})


if __name__ == "__main__":
    import asyncio

    from pyeed.database import Database

    async def main() -> None:
        db = Database("bolt://129.69.129.132:7687", "neo4j", "12345678")

        await embed_proteins(
            db,
            model_name="facebook/esm2_t33_650M_UR50D",
            gpu_batch_size=26,
            stream_chunk_size=104,  # Use reasonable chunk size to avoid memory thrashing
            db_batch_size=208,  # Batch DB writes for efficiency
        )

        await db.close()

    asyncio.run(main())
