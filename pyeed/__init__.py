from .db.neo4j import Database
from .environment import IN_IPYTHON, IN_NOTEBOOK
from .ingest import (
    enrich_molecules,
    enrich_reactions,
    fetch_proteins_by_ids,
    fetch_proteins_by_interpro,
    ingest_full_pipeline,
)
from .logging_setup import setup_logging

setup_logging()

__all__ = [
    "IN_IPYTHON",
    "IN_NOTEBOOK",
    "Database",
    "enrich_molecules",
    "enrich_reactions",
    "fetch_proteins_by_ids",
    "fetch_proteins_by_interpro",
    "ingest_full_pipeline",
]
