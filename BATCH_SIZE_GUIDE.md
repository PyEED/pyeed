# Batch Size Configuration Guide

## Single Source of Truth ✅

All batch sizes are now clearly defined in `embed_proteins()` with distinct purposes:

```python
await embed_proteins(
    db,
    gpu_batch_size=8,        # GPU memory constraint
    stream_chunk_size=1000,  # Progress update frequency
    db_batch_size=100,       # Database transaction size
)
```

## Three Batch Size Concepts

### 1. `gpu_batch_size` (default: 8)
**Purpose**: Controls how many sequences are processed per GPU batch  
**Located**: `esm2_async.py` → `embed_batch()` → `_compute_embeddings()`  
**Constraint**: GPU memory  

**Tune based on:**
- GPU VRAM (higher = more memory needed)
- Sequence length (longer sequences = more memory)
- Model size (ESM-2 650M vs 3B)

**Examples:**
- Small GPU (8GB): `gpu_batch_size=4`
- Medium GPU (16GB): `gpu_batch_size=8`
- Large GPU (40GB+): `gpu_batch_size=16-32`

### 2. `stream_chunk_size` (default: 1000)
**Purpose**: Controls how often embeddings are yielded and progress updates  
**Located**: `esm2_async.py` → `embed_stream()`  
**Constraint**: Progress bar refresh frequency  

**Tune based on:**
- How often you want progress updates (lower = more frequent)
- Memory overhead (higher = fewer yields but more memory)
- Responsiveness vs throughput trade-off

**Examples:**
- Fast updates: `stream_chunk_size=100` (updates every ~10 seconds)
- Balanced: `stream_chunk_size=1000` (updates every ~1-2 minutes)
- Throughput: `stream_chunk_size=5000` (updates every ~5-10 minutes)

### 3. `db_batch_size` (default: 100)
**Purpose**: Controls how many proteins are accumulated before database write  
**Located**: `embed_proteins.py` → `consumer()`  
**Constraint**: Database transaction overhead  

**Tune based on:**
- Database performance (larger batches = fewer transactions)
- Memory for queue (larger = more memory)
- Failure recovery (smaller = less work lost on error)

**Examples:**
- Conservative: `db_batch_size=50`
- Balanced: `db_batch_size=100`
- Aggressive: `db_batch_size=500`

## Data Flow

```
Input: 200,000 sequences
   ↓
┌─────────────────────────────────────────────────┐
│ embed_stream() chunks by stream_chunk_size      │
│ Chunks: 200,000 / 1000 = 200 chunks             │
└─────────────────────────────────────────────────┘
   ↓ (yields every 1000)
┌─────────────────────────────────────────────────┐
│ embed_batch() processes gpu_batch_size at a time│
│ Per chunk: 1000 / 8 = 125 GPU batches           │
│ Multi-GPU: Distributed round-robin              │
└─────────────────────────────────────────────────┘
   ↓ (streams to queue)
┌─────────────────────────────────────────────────┐
│ Queue: maxsize = db_batch_size * 2              │
│ Buffers between producer and consumer           │
└─────────────────────────────────────────────────┘
   ↓ (consumed)
┌─────────────────────────────────────────────────┐
│ consumer() accumulates db_batch_size proteins   │
│ Then: db.attach() writes batch to Neo4j         │
└─────────────────────────────────────────────────┘
```

## Progress Updates

Progress bars update at **`stream_chunk_size`** intervals:

```
[ESM2] Completed chunk 0-1000 (1000/203355)
Progress: Compute Embeddings ████░░░░ 1000/203355

[ESM2] Completed chunk 1000-2000 (2000/203355)
Progress: Compute Embeddings ████████░ 2000/203355
```

**Why progress was stuck at 0:**
- Progress updates happened in the async loop
- But `embed_stream()` used synchronous `for` loop
- Fixed by ensuring proper async yielding

## Recommended Configurations

### Small Dataset (<10K sequences)
```python
gpu_batch_size=8,
stream_chunk_size=100,   # Frequent updates
db_batch_size=50,
```

### Medium Dataset (10K-100K sequences)
```python
gpu_batch_size=16,
stream_chunk_size=1000,  # Balanced
db_batch_size=100,
```

### Large Dataset (>100K sequences)
```python
gpu_batch_size=16,
stream_chunk_size=5000,  # Less frequent updates
db_batch_size=500,       # Larger DB batches
```

### Multi-GPU Setup (4x GPUs)
```python
gpu_batch_size=16,       # Per GPU
stream_chunk_size=1000,  # 1000 / 4 GPUs = 250 per GPU
db_batch_size=200,
```

## Memory Estimates

For 200K sequences with ESM-2 650M:

**GPU Memory (per device):**
- Model: ~2.5 GB
- Batch (8 seqs): ~1-2 GB
- Peak: ~4-5 GB per GPU

**System Memory:**
- Queue: `db_batch_size * 2 * embedding_size`
  - 100 * 2 * 1280 * 4 bytes ≈ 1 MB
- Chunk buffer: `stream_chunk_size * embedding_size`
  - 1000 * 1280 * 4 bytes ≈ 5 MB
- Total: <100 MB overhead

## Troubleshooting

### Progress bar stuck at 0
- **Cause**: Async loop not yielding properly
- **Fix**: Applied in this commit (proper async yielding)

### OOM errors
- **Solution**: Reduce `gpu_batch_size`
- Start at 4, increase gradually

### Slow progress updates
- **Solution**: Reduce `stream_chunk_size`
- More frequent updates = more overhead

### Database transaction errors
- **Solution**: Reduce `db_batch_size`
- Smaller batches = more transactions but more reliable

