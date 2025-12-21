"""Protein embedding module.

Provides high-performance ESM2 protein embeddings with multi-device support.

Core components:
    ESM2Processor: Multi-device ESM2 embedder with length-sorted batching.
    run_embedding_job: End-to-end pipeline from Neo4j to Milvus.

Example:
    >>> from pyeed.embed import ESM2Processor, run_embedding_job
    >>> 
    >>> # Standalone embedding
    >>> async with ESM2Processor(devices=[0], dtype="float16") as processor:
    ...     results = await processor.work([("P12345", "MVLSPADKTN...")])
    >>> 
    >>> # Full pipeline
    >>> await run_embedding_job(
    ...     neo4j_driver=driver,
    ...     milvus_client=client,
    ...     collection_name="proteins",
    ... )
"""

from .embed_from_db import run_embedding_job
from .embedder import ESM2Processor

__all__ = [
    "ESM2Processor",
    "run_embedding_job",
]

