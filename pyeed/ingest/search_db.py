from ..db.milvus import VectorDB

def ingest_fasta(vector_db: VectorDB, path: str) -> None:
    