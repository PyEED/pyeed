from .database import Database
from .logging_config import setup_logging
from .main import ingest_uniprot

setup_logging()

__all__ = ["Database", "ingest_uniprot"]
