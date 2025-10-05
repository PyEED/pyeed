"""
Async protein embedding pipeline with streaming database writes.

Processes large protein datasets efficiently by:
- Finding proteins without embeddings
- Computing embeddings in batches
- Streaming results to database as they complete
"""

import asyncio

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
from pyeed.embeddings.pooling import mean_pooling
from pyeed.model import Protein
from pyeed.model.embedding import Embedding


async def get_proteins_without_embeddings(
    db: Database,
    model_name: str,
    pooling_method: str = "mean_pooling",
    layer_index: int = -1,
) -> list[tuple[str, str]]:
    """
    Query proteins that don't have embeddings for specified model/pooling/layer.

    Returns:
        List of (accession_id, sequence) tuples
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
    RETURN p.accession_id AS accession_id, p.sequence AS sequence
    """

    print(
        f"[Embed] Querying proteins without embeddings (model={model_name}, pooling={pooling_method}, layer={layer_index})..."
    )

    results = []
    async for record in db.async_query_iter(
        query,
        model_name=model_name_clean,
        pooling_method=pooling_clean,
        layer_index=layer_index,
    ):
        results.append((record["accession_id"], record["sequence"]))

    print(f"[Embed] Found {len(results)} proteins without embeddings")
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
        print("[Embed] No proteins to embed")
        return

    accessions = [acc for acc, _ in proteins_todo]
    sequences = [seq for _, seq in proteins_todo]

    print(f"[Embed] Starting embedding pipeline for {len(proteins_todo)} proteins")
    print(
        f"[Embed] GPU batch: {gpu_batch_size}, Stream chunks: {stream_chunk_size}, DB batch: {db_batch_size}"
    )

    # Create queue for streaming results: (accession_id, Embedding)
    queue: asyncio.Queue[tuple[str, Embedding] | None] = asyncio.Queue(maxsize=db_batch_size * 2)

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
                # Update progress for entire chunk at once (more efficient)
                chunk_size = len(embeddings_chunk)
                for i, emb in enumerate(embeddings_chunk):
                    acc = accessions[chunk_start + i]
                    await queue.put((acc, emb))

                # Update progress once per chunk instead of per embedding
                progress.update(task_id, advance=chunk_size)
                progress.refresh()  # Force refresh
                chunk_start += chunk_size

        # Signal done
        await queue.put(None)

    async def consumer(task_id: TaskID) -> None:
        """Attach embeddings to existing proteins in database."""
        # Build batches: accession_id -> [Embedding, ...]
        batch: dict[str, list[Embedding]] = {}

        while True:
            item = await queue.get()
            if item is None:
                # Final flush
                if batch:
                    await db.attach(Protein, batch)
                    progress.update(task_id, advance=sum(len(v) for v in batch.values()))
                break

            # Unpack (accession, embedding) tuple
            accession_id, embedding = item

            # Group embeddings by protein accession
            if accession_id not in batch:
                batch[accession_id] = []
            batch[accession_id].append(embedding)

            # Save batch when it reaches target size
            if len(batch) >= db_batch_size:
                await db.attach(Protein, batch)
                count = sum(len(v) for v in batch.values())
                progress.update(task_id, advance=count)
                progress.refresh()  # Force refresh
                batch.clear()

    with progress:
        embed_task = progress.add_task("[cyan]Compute Embeddings", total=len(proteins_todo))
        save_task = progress.add_task("[green]Save to DB", total=len(proteins_todo))

        await asyncio.gather(
            producer(embed_task),
            consumer(save_task),
        )

    print(f"[Embed] Complete! Embedded {len(proteins_todo)} proteins")


if __name__ == "__main__":
    import asyncio

    from pyeed.database import Database

    async def main() -> None:
        db = Database("bolt://129.69.129.132:7687", "neo4j", "12345678")

        await embed_proteins(
            db,
            model_name="facebook/esm2_t33_650M_UR50D",
            gpu_batch_size=8,
            stream_chunk_size=8,
            db_batch_size=8,
        )

        await db.close()

    asyncio.run(main())
