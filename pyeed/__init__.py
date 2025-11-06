from .db.neo4j import Database
from .environment import IN_IPYTHON, IN_NOTEBOOK
from .ingest import (
    Ingester,
    enrich_molecules,
    enrich_reactions,
    fetch_proteins_by_ids,
    fetch_proteins_by_interpro,
    ingest_full_pipeline,
)
from .logging_setup import setup_logging

# Backward compatibility (deprecated)
from .main import enrich_rhea, ingest_interpro, ingest_uniprot

setup_logging()

__all__ = [
    "IN_IPYTHON",
    "IN_NOTEBOOK",
    "Database",
    "Ingester",
    "enrich_molecules",
    "enrich_reactions",
    "enrich_rhea",
    "fetch_proteins_by_ids",
    "fetch_proteins_by_interpro",
    "ingest_full_pipeline",
    "ingest_interpro",
    "ingest_uniprot",
]
