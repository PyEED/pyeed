from __future__ import annotations

import asyncio
from collections.abc import Callable
from typing import Literal

import pandas as pd
from loguru import logger
from rich.console import Group
from rich.live import Live
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from pyeed.db.milvus import VectorDB
from pyeed.embed.esm2 import ESM2Embedder
from pyeed.embed.pooling import PoolingLike, mean_pooling

from .db.neo4j import GraphDB
from .ingest.core.pipeline import Pipeline
from .ingest.sources.fasta import build_header_index
from .ingest.stages.embed import EmbeddingStage
from .ingest.stages.fasta import FASTAReaderStage
from .ingest.stages.sink import Neo4jUpsertStage
from .utils.progress import CONSOLE


def ingest_fasta(
    fasta_path: str,
    chunk_size: int = 10,
    header_fn: Callable[[str], str] | None = None,
    graph_db_uri: str | None = None,
    graph_db_user: str | None = None,
    graph_db_password: str | None = None,
    vector_db_uri: str | None = None,
    vector_db_token: str | None = None,
    vector_db_collection: str = "pyeed",
    huggingface_token: str | None = None,
    huggingface_model_name: str = "facebook/esm2_t33_650M_UR50D",
    huggingface_model_max_seq_length: int = 1024,
    model_dtype: Literal["float16", "float32"] = "float32",
    return_dtype: Literal["float16", "float32"] = "float32",
    n_gpus: int = -1,
    pooling_methods: PoolingLike = [mean_pooling],
    embedding_chunk_size: int = 1000,
    embedding_batch_size: int = 32,
    enable_embedding: bool = True,
) -> None:
    """Ingest a FASTA file into the pipeline.

    Args:
        fasta_path: Path to FASTA file
        chunk_size: Number of sequences per chunk for FASTA reading
        header_fn: Function to extract protein_id from header
            Signature: (header: str) -> str
            If None, uses entire header as id of protein
        graph_db_uri: GraphDB connection URI
        graph_db_user: GraphDB username
        graph_db_password: GraphDB password
        vector_db_uri: VectorDB (Milvus) connection URI
        vector_db_token: VectorDB authentication token
        vector_db_collection: Collection name in VectorDB
        huggingface_token: HuggingFace API token (optional, for gated models)
        huggingface_model_name: ESM2 model name
        huggingface_model_max_seq_length: Maximum sequence length
        model_dtype: Model computation dtype
        return_dtype: Embedding return dtype
        n_gpus: Number of GPUs to use (-1 for all available)
        pooling_methods: Pooling methods for embeddings
        embedding_chunk_size: Records to accumulate before embedding (default: 1000)
        embedding_batch_size: GPU batch size for embedding (default: 32)
        enable_embedding: Whether to enable embedding stage (default: True)

    Returns:
        None
    """
    bar = Progress(
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        TimeElapsedColumn(),
        transient=False,
        # console=CONSOLE,
    )
    spin = Progress(
        SpinnerColumn(style="cyan", spinner_name="dots", finished_text="✓"),
        TextColumn("{task.description}"),
        transient=False,
        # console=CONSOLE,
    )

    # Tasks
    graph_db_connect_task = spin.add_task("Connecting to GraphDB", total=1)
    vector_db_connect_task = spin.add_task("Connecting to VectorDB", total=1)
    embedder_init_task = spin.add_task("Initializing ESM2Embedder", total=1)
    read_task = bar.add_task("Read FASTA", total=None)
    embed_task = bar.add_task("Embed Sequences", total=None)
    write_task = bar.add_task("Add Sequences to GraphDB", total=None)

    with Live(Group(spin, bar), console=CONSOLE, refresh_per_second=8):
        # Connect to GraphDB
        graph_db = GraphDB(
            uri=graph_db_uri,
            user=graph_db_user,
            password=graph_db_password,
        )
        logger.info(f"Connected to GraphDB: {graph_db.uri}")
        spin.update(graph_db_connect_task, completed=1)
        spin.update(graph_db_connect_task, description=f"Connected to GraphDB: {graph_db.uri}")

        # Connect to VectorDB
        vector_db = VectorDB(
            uri=vector_db_uri,
            token=vector_db_token,
        )
        spin.update(vector_db_connect_task, completed=1)
        spin.update(vector_db_connect_task, description=f"Connected to VectorDB: {vector_db.uri}")
        logger.info(f"Connected to VectorDB: {vector_db.uri}")

        # Initialize ESM2Embedder (only if embedding is enabled)
        embedder = None
        if enable_embedding:
            embedder = ESM2Embedder(
                model_name=huggingface_model_name,
                model_dtype=model_dtype,
                return_dtype=return_dtype,
                n_gpus=n_gpus,
                pooling_methods=pooling_methods,
                max_length=huggingface_model_max_seq_length,
            )
            asyncio.run(embedder.initialize())
            spin.update(embedder_init_task, completed=1)
            spin.update(
                embedder_init_task,
                description=f"Initialized {embedder.model_name} on {len(embedder.devices)} devices",
            )
            logger.info("Initialized ESM2Embedder")
        else:
            spin.update(embedder_init_task, completed=1)
            spin.update(embedder_init_task, description="Embedding disabled")
            logger.info("Embedding stage disabled")

        # Build header index from FASTA file
        offsets = build_header_index(fasta_path)
        total_sequences = len(offsets)
        bar.update(read_task, description="Read FASTA", total=total_sequences)
        bar.update(embed_task, description="Embed Sequences", total=total_sequences)
        bar.update(write_task, description="Add Sequences to GraphDB", total=total_sequences)

        # Build pipeline with optional embedding stage
        pipeline = Pipeline(progress=bar)

        if enable_embedding and embedder is not None:
            # Three-stage pipeline: Read -> Embed -> Write
            # Queue sizes based on chunk_size with buffer for backpressure
            queue_size = max(embedding_chunk_size + 100, 1000)
            read_queue = pipeline.add_queue("read", maxsize=queue_size)
            embed_queue = pipeline.add_queue("embed", maxsize=queue_size)

            pipeline.add_stage(
                stage=FASTAReaderStage(
                    fasta_path=fasta_path,
                    offsets=offsets,
                    chunk_size=chunk_size,
                    header_extractor=header_fn,
                ),
                input_queues=[],
                output_queues=[read_queue],
                task_id=read_task,
            )

            pipeline.add_stage(
                stage=EmbeddingStage(
                    embedder=embedder,
                    chunk_size=embedding_chunk_size,
                    batch_size=embedding_batch_size,
                ),
                input_queues=[read_queue],
                output_queues=[embed_queue],
                task_id=embed_task,
            )

            pipeline.add_stage(
                stage=Neo4jUpsertStage(
                    db=graph_db,
                    batch_size=100,
                ),
                input_queues=[embed_queue],
                output_queues=[],
                task_id=write_task,
            )
        else:
            # Two-stage pipeline: Read -> Write (no embedding)
            node_queue = pipeline.add_queue("nodes", maxsize=10)

            pipeline.add_stage(
                stage=FASTAReaderStage(
                    fasta_path=fasta_path,
                    offsets=offsets,
                    chunk_size=chunk_size,
                    header_extractor=header_fn,
                ),
                input_queues=[],
                output_queues=[node_queue],
                task_id=read_task,
            )

            pipeline.add_stage(
                stage=Neo4jUpsertStage(
                    db=graph_db,
                    batch_size=100,
                ),
                input_queues=[node_queue],
                output_queues=[],
                task_id=write_task,
            )

        asyncio.run(pipeline.run())


def ingest_csv(
    path: str | pd.DataFrame,
    graph_db_uri: str | None = None,
    graph_db_user: str | None = None,
    graph_db_password: str | None = None,
    vector_db_uri: str | None = None,
    vector_db_token: str | None = None,
    vector_db_collection: str = "pyeed",
) -> None:
    """Ingest a CSV file into the pipeline."""
    pass


def ingest_uniprot(
    ids: list[str],
    graph_db_uri: str | None = None,
    graph_db_user: str | None = None,
    graph_db_password: str | None = None,
    vector_db_uri: str | None = None,
    vector_db_token: str | None = None,
    vector_db_collection: str = "pyeed",
) -> None:
    """Ingest a list of UniProt IDs into the pipeline."""
    pass


def ingest_interpro(
    id: str,
    graph_db_uri: str | None = None,
    graph_db_user: str | None = None,
    graph_db_password: str | None = None,
    vector_db_uri: str | None = None,
    vector_db_token: str | None = None,
    vector_db_collection: str = "pyeed",
) -> None:
    """Ingest an InterPro ID into the pipeline."""
    pass


if __name__ == "__main__":
    ingest_fasta(fasta_path="/home/mha/projects/proteingraph/downloads/1000seq.fasta", n_gpus=2)
