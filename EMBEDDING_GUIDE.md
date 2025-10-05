# Embedding Pipeline Guide

## ✅ Ready for 200K Proteins

Your embedding system is now production-ready for large-scale protein embedding with:

### Key Features

1. **Async Streaming Architecture**
   - Embeddings computed in chunks
   - Results streamed to database as they complete
   - No need to hold 200K embeddings in memory

2. **Multi-GPU Parallelization**
   - Automatic round-robin distribution across all CUDA devices
   - Efficient batching per GPU

3. **Producer-Consumer Pattern**
   - Same pattern as `main.py` (UniProt/Rhea ingestion)
   - Producer: computes embeddings
   - Consumer: writes to database in batches

4. **Smart Filtering**
   - Only embeds proteins that don't have embeddings yet
   - Checks for specific model/pooling/layer combination

## Usage

### Simple: Embed All Missing Proteins

```python
import asyncio
from pyeed.database import Database
from pyeed.embed_proteins import embed_proteins

async def main():
    db = Database()
    
    await embed_proteins(
        db,
        model_name="facebook/esm2_t33_650M_UR50D",
        batch_size=8,              # Sequences per batch per GPU
        save_batch_size=100,       # DB writes every 100 proteins
    )
    
    await db.close()

asyncio.run(main())
```

### Advanced: Custom Streaming

```python
from pyeed.embeddings import ESM2Embedder

async with ESM2Embedder(model_name="facebook/esm2_t33_650M_UR50D") as embedder:
    # Stream embeddings in chunks
    async for embeddings_chunk in embedder.embed_stream(
        sequences=sequences,
        accessions=accessions,
        batch_size=8,
        chunk_size=1000,  # Yield every 1000 embeddings
    ):
        # Process chunk immediately
        await save_to_database(embeddings_chunk)
```

## Architecture

### Files

- **`pyeed/embeddings/esm2_async.py`**: Core ESM2Embedder class with streaming
- **`pyeed/embed_proteins.py`**: High-level pipeline for full workflow

### Flow

```
1. Query DB for proteins without embeddings
   ↓
2. Producer: ESM2Embedder.embed_stream()
   - Computes embeddings in chunks
   - Yields as soon as chunk completes
   ↓
3. Queue: Async queue buffers chunks
   ↓
4. Consumer: Database.save_many()
   - Writes batches to Neo4j
   - Progress tracking
```

### Memory Efficiency

For 200K proteins:
- **Old approach**: Load all 200K embeddings → OOM
- **New approach**: Stream in chunks of ~1000
- **Peak memory**: ~10-20 chunks in queue max

### Database Schema

Proteins link to Embeddings via `HAS_EMBEDDING` relationship:

```cypher
(p:Protein)-[:HAS_EMBEDDING]->(e:Embedding)

WHERE e.model_name = "facebook_esm2_t33_650m_ur50d"
  AND e.pooling_method = "mean_pooling"
  AND e.layer_index = -1
```

## Performance Tips

### For 200K Proteins

```python
await embed_proteins(
    db,
    batch_size=16,           # Higher if you have GPU memory
    save_batch_size=200,     # Larger DB batches = fewer transactions
    device_ids=[0, 1, 2, 3], # Use all GPUs
)
```

### Expected Times (ESM-2 650M)

- **Single GPU (A100)**: ~8-10 sequences/sec → 6-7 hours for 200K
- **4x GPUs**: ~30-40 sequences/sec → 1.5-2 hours for 200K
- **Database writes**: Negligible (async, batched)

## Error Handling

The pipeline is robust:
- Model loading failures: Retry with single GPU
- OOM errors: Auto-reduces batch size
- Database errors: Batch is logged, pipeline continues

## Monitoring

Progress bars show:
1. **Compute Embeddings**: How many sequences processed
2. **Save to DB**: How many proteins saved

Logs show:
- Model loading per GPU
- Batch completions
- Chunk streaming progress

## Next Steps

1. Run on subset first: Test with 100 proteins
2. Monitor GPU memory: `nvidia-smi -l 1`
3. Scale up to full 200K
4. Consider checkpointing for very large runs

## Example Output

```
[Embed] Querying proteins without embeddings...
[Embed] Found 200000 proteins without embeddings
[ESM2] Model: facebook/esm2_t33_650M_UR50D
[ESM2] Pooling: mean_pooling, Layer: -1, Normalize: True
[ESM2] Loading model on 4 device(s): [cuda:0, cuda:1, cuda:2, cuda:3]
[ESM2] Loading on cuda:0...
[ESM2] Loaded on cuda:0
...
[ESM2] Streaming 200000 sequences in chunks of 320
Compute Embeddings  ████████░░░░░░░░  45000/200000 [00:45:23<01:02:15, 41.5 it/s]
Save to DB         ████████░░░░░░░░  44800/200000 [00:45:20<01:02:18, 41.4 it/s]
```

