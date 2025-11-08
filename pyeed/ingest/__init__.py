"""Data ingestion and enrichment API for pyeed.

This module provides both functional and fluent interfaces for ingesting
protein data from UniProt and enriching it with Rhea reactions, ChEBI molecules,
and protein embeddings.

Functional API:
    Use individual async functions for fine-grained control:
    - fetch_proteins_by_ids()
    - fetch_proteins_by_interpro()
    - enrich_reactions()
    - enrich_molecules()
    - embed_proteins()
    - ingest_full_pipeline()

Fluent API:
    Use the Ingester class for chainable configuration:
    - Ingester(db).from_interpro("IPR002133").with_embeddings().run()

Embedding Sinks:
    - EmbeddingSink (Protocol)
"""

from .functions import (
    enrich_molecules,
    enrich_reactions,
    fetch_proteins_by_ids,
    fetch_proteins_by_interpro,
    ingest_full_pipeline,
)
from .pipeline import EmbeddingPipeline

__all__ = [
    "EmbeddingPipeline",
    "enrich_molecules",
    "enrich_reactions",
    "fetch_proteins_by_ids",
    "fetch_proteins_by_interpro",
    "ingest_full_pipeline",
]
