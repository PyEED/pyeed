# Unified Pipeline Orchestrator Refactor

## Overview

Create a flexible, stage-based async pipeline that can:

- Ingest from 3 sources: FASTA, UniProt IDs/InterPro, CSV
- Optionally enrich with Rhea reactions and ChEBI molecules
- Optionally embed sequences with ESM2
- Write metadata to Neo4j and vectors to Milvus
- Skip existing entries (check both databases)
- Maintain backpressure via bounded queues

## Architecture

### Stage-Based Design (Producer-Consumer Pattern)

```
Source Stage → Enrichment Stage → Embedding Stage → Sink Stage(s)
     |               |                   |               |
     └─ queue ──────→└─ queue ──────────→└─ queue ──────→
```

Each stage is an async worker with bounded input/output queues. The embedder controls overall throughput via backpressure.

### Core Components

#### 1. Data Model: `PipelineRecord`

Create `pyeed/ingest/pipeline_types.py`:

- `PipelineRecord`: Container for data flowing through pipeline
  - `sequence_id: str` (primary key)
  - `sequence: str`
  - `metadata: dict[str, Any]` (flexible for CSV/UniProt extra data)
  - `protein: Protein | None` (Pydantic model, if enriched)
  - `embeddings: dict[str, np.ndarray] | None` (if embedded)
  - `skip_neo4j: bool`, `skip_milvus: bool` (flags for existence checks)

#### 2. Source Adapters

Create `pyeed/ingest/sources/`:

**`base.py`**:

- `Protocol`: `PipelineSource` with `async def stream() -> AsyncIterator[PipelineRecord]`

**`fasta_source.py`**:

- `FASTASource(PipelineSource)`: Reuse existing `read_fasta_chunks_async`
- Emit `PipelineRecord(sequence_id, sequence, metadata={})`

**`uniprot_source.py`**:

- `UniProtSource(PipelineSource)`: Use existing `UniProtAdapter`
- Support both accession list and InterPro ID (calls SPARQL first)
- Emit `PipelineRecord` with `protein` field populated
- Fetch in chunks with configurable concurrency

**`csv_source.py`**:

- `CSVSource(PipelineSource)`: Read CSV with required columns: `sequence_id`, `sequence`
- Optional columns validated against Pydantic models if present
- Extra columns go into `metadata` dict freestyle
- Support async chunk reading for large files

#### 3. Transform Stages

Create `pyeed/ingest/stages/`:

**`base.py`**:

- `Protocol`: `PipelineStage` with `async def process(record: PipelineRecord) -> PipelineRecord`

**`existence_check_stage.py`**:

- `ExistenceCheckStage`: Check Neo4j (by sequence_id) and Milvus (by collection)
- Set `record.skip_neo4j` and `record.skip_milvus` flags
- Runs in thread pool to avoid blocking
- Batched checks for efficiency (check 100s at once)

**`enrichment_stage.py`**:

- `EnrichmentStage`: If `record.protein` exists and has `rhea_id`, fetch reactions/molecules
- Uses existing `RheaClient` and `ChebiClient`
- Controlled concurrency via semaphore
- Updates `record.protein` with enriched data

**`embedding_stage.py`**:

- `EmbeddingStage`: Wrapper around existing `ESM2Embedder`
- Accumulates records into batches (sort by length)
- Calls `embedder.embed_stream()` with backpressure
- Populates `record.embeddings` dict
- Handles long sequences (zero embeddings, like current implementation)

#### 4. Sink Stages

**`neo4j_sink.py`**:

- `Neo4jSink`: Batch writes to Neo4j
- Extract sequence_id, sequence, seq_length + metadata
- Use existing `db.save_many()` for Protein nodes
- Respect `record.skip_neo4j` flag

**`milvus_sink.py`**:

- `MilvusSink`: Batch writes to Milvus
- Use existing `vector_db.insert_async()`
- Respect `record.skip_milvus` flag

#### 5. Orchestrator

Create `pyeed/ingest/orchestrator.py`:

**`PipelineOrchestrator`**:

- Constructor takes: source, stages list, sinks list, queue sizes
- `async def run() -> PipelineStats`: Main entry point
- Spawns async workers for each stage with queues between them
- Uses `SENTINEL` pattern for graceful shutdown
- Progress tracking with `rich.progress` (unified view of all stages)
- Returns statistics (records processed, skipped, errors)

**Configuration via params** (not separate config class for simplicity):

```python
async def run_pipeline(
    # Source selection (exactly one required)
    fasta_path: str | None = None,
    uniprot_ids: list[str] | None = None,
    interpro_id: str | None = None,
    csv_path: str | None = None,
    # Databases
    neo4j_db: Database,
    vector_db: VectorDB | None = None,  # None = skip embedding
    milvus_collection: str | None = None,
    # Embedder (optional)
    embedder: ESM2Embedder | None = None,
    # Enrichment flags
    enrich_reactions: bool = False,
    enrich_molecules: bool = False,
    # Performance tuning
    chunk_size: int = 3200,
    batch_size: int = 32,
    write_batch_size: int = 10,
    queue_size: int = 256,
    # Behavior
    skip_existing: bool = True,
    header_extractor: Callable[[str], str] | None = None,  # FASTA only
) -> PipelineStats:
    ...
```

## File Structure

```
pyeed/ingest/
├── pipeline.py              # Legacy - keep for now, deprecate later
├── pipeline_types.py        # NEW: PipelineRecord, PipelineStats
├── orchestrator.py          # NEW: PipelineOrchestrator, run_pipeline()
├── sources/
│   ├── __init__.py
│   ├── base.py             # NEW: PipelineSource protocol
│   ├── fasta_source.py     # NEW
│   ├── uniprot_source.py   # NEW
│   └── csv_source.py       # NEW
├── stages/
│   ├── __init__.py
│   ├── base.py             # NEW: PipelineStage protocol
│   ├── existence_check_stage.py  # NEW
│   ├── enrichment_stage.py       # NEW
│   └── embedding_stage.py        # NEW
├── sinks/
│   ├── __init__.py
│   ├── base.py             # NEW: PipelineSink protocol
│   ├── neo4j_sink.py       # NEW
│   └── milvus_sink.py      # NEW
├── functions.py            # Legacy - migrate to orchestrator gradually
├── uniprot.py              # Keep as adapter
├── rhea.py                 # Keep as client
├── chebi.py                # Keep as client
└── progress.py             # Keep as utility
```

## Implementation Strategy

1. Create type definitions and protocols (foundation)
2. Implement source adapters (reuse existing code)
3. Implement stages (reuse embedder, clients)
4. Implement sinks (reuse DB clients)
5. Implement orchestrator (wire everything together)
6. Add comprehensive logging and error handling
7. Create usage examples for each source type
8. Update module `__init__.py` to expose new API

## Key Design Principles

- **Separation of concerns**: Each stage has single responsibility
- **Dependency injection**: Pass clients/embedder to orchestrator
- **Protocol-based**: Use `Protocol` for extensibility
- **Backpressure**: Bounded queues propagate upstream
- **Fail-fast validation**: Check inputs early
- **Observability**: Structured logging + progress bars
- **Testing**: Each stage testable independently
- **No global state**: All state in orchestrator instance

## Migration Path

- Keep existing `pipeline.py` and `functions.py` working
- New code in separate modules
- Users can adopt incrementally
- Mark old code as deprecated in docstrings
- Remove old code in future version