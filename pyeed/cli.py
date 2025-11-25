"""CLI tool for pyeed pipelines using Typer."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal

import typer

from .pipelines import ingest_fasta, ingest_interpro, ingest_uniprot

app = typer.Typer(
    name="pyeed",
    help="Pyeed: Toolkit to create, annotate, and analyze sequence data",
    add_completion=True,
)


@app.command(name="ingest-fasta")
def ingest_fasta_cmd(
    fasta_path: Annotated[
        Path,
        typer.Argument(
            help="Path to FASTA file to ingest",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
        ),
    ],
    chunk_size: Annotated[
        int,
        typer.Option(
            "--chunk-size",
            help="Number of sequences per chunk for FASTA reading",
        ),
    ] = 512,
    graph_db_uri: Annotated[
        str | None,
        typer.Option(
            "--graph-db-uri",
            envvar="NEO4J_URI",
            help="Neo4j connection URI",
        ),
    ] = None,
    graph_db_user: Annotated[
        str | None,
        typer.Option(
            "--graph-db-user",
            envvar="NEO4J_USER",
            help="Neo4j username",
        ),
    ] = None,
    graph_db_password: Annotated[
        str | None,
        typer.Option(
            "--graph-db-password",
            envvar="NEO4J_PASSWORD",
            help="Neo4j password",
        ),
    ] = None,
    vector_db_uri: Annotated[
        str | None,
        typer.Option(
            "--vector-db-uri",
            envvar="MILVUS_URI",
            help="Milvus connection URI",
        ),
    ] = None,
    vector_db_token: Annotated[
        str | None,
        typer.Option(
            "--vector-db-token",
            envvar="MILVUS_TOKEN",
            help="Milvus authentication token",
        ),
    ] = None,
    vector_db_collection: Annotated[
        str,
        typer.Option(
            "--vector-db-collection",
            help="Collection name in VectorDB",
        ),
    ] = "pyeed",
    huggingface_token: Annotated[
        str | None,
        typer.Option(
            "--huggingface-token",
            envvar="HUGGINGFACE_TOKEN",
            help="HuggingFace API token (for gated models)",
        ),
    ] = None,
    model_name: Annotated[
        str,
        typer.Option(
            "--model-name",
            help="ESM2 model name from HuggingFace",
        ),
    ] = "facebook/esm2_t33_650M_UR50D",
    n_gpus: Annotated[
        int,
        typer.Option(
            "--n-gpus",
            help="Number of GPUs to use (-1 for all available)",
        ),
    ] = -1,
    model_dtype: Annotated[
        Literal["float16", "float32"],
        typer.Option(
            "--model-dtype",
            help="Model computation dtype",
        ),
    ] = "float32",
    return_dtype: Annotated[
        Literal["float16", "float32"],
        typer.Option(
            "--return-dtype",
            help="Embedding return dtype",
        ),
    ] = "float32",
    embedding_batch_size: Annotated[
        int,
        typer.Option(
            "--embedding-batch-size",
            help="GPU batch size for embedding",
        ),
    ] = 32,
    embedding_chunk_size: Annotated[
        int,
        typer.Option(
            "--embedding-chunk-size",
            help="Records to accumulate before embedding",
        ),
    ] = 1000,
    enrichment_batch_size: Annotated[
        int,
        typer.Option(
            "--enrichment-batch-size",
            help="Records to accumulate before enriching",
        ),
    ] = 50,
    enrichment_max_concurrent: Annotated[
        int,
        typer.Option(
            "--enrichment-max-concurrent",
            help="Maximum concurrent API requests for enrichment",
        ),
    ] = 50,
) -> None:
    """Ingest a FASTA file into the pipeline.

    Reads sequences from a FASTA file, embeds them using ESM2, and stores
    them in Neo4j (graph database) and Milvus (vector database).
    """
    ingest_fasta(
        fasta_path=str(fasta_path),
        chunk_size=chunk_size,
        header_fn=None,
        taxon_fn=None,
        graph_db_uri=graph_db_uri,
        graph_db_user=graph_db_user,
        graph_db_password=graph_db_password,
        vector_db_uri=vector_db_uri,
        vector_db_token=vector_db_token,
        vector_db_collection=vector_db_collection,
        huggingface_token=huggingface_token,
        huggingface_model_name=model_name,
        huggingface_model_max_seq_length=1024,
        model_dtype=model_dtype,
        return_dtype=return_dtype,
        n_gpus=n_gpus,
        pooling_methods=None,
        embedding_chunk_size=embedding_chunk_size,
        embedding_batch_size=embedding_batch_size,
        enrichment_batch_size=enrichment_batch_size,
        enrichment_max_concurrent=enrichment_max_concurrent,
    )


@app.command(name="ingest-uniprot")
def ingest_uniprot_cmd(
    ids: Annotated[
        list[str],
        typer.Argument(
            help="List of UniProt accession IDs (e.g., P12345 Q9Y6X9)",
        ),
    ],
    uniprot_chunk_size: Annotated[
        int,
        typer.Option(
            "--uniprot-chunk-size",
            help="Accessions per API batch",
        ),
    ] = 50,
    uniprot_page_size: Annotated[
        int,
        typer.Option(
            "--uniprot-page-size",
            help="Results per page",
        ),
    ] = 50,
    graph_db_uri: Annotated[
        str | None,
        typer.Option(
            "--graph-db-uri",
            envvar="NEO4J_URI",
            help="Neo4j connection URI",
        ),
    ] = None,
    graph_db_user: Annotated[
        str | None,
        typer.Option(
            "--graph-db-user",
            envvar="NEO4J_USER",
            help="Neo4j username",
        ),
    ] = None,
    graph_db_password: Annotated[
        str | None,
        typer.Option(
            "--graph-db-password",
            envvar="NEO4J_PASSWORD",
            help="Neo4j password",
        ),
    ] = None,
    vector_db_uri: Annotated[
        str | None,
        typer.Option(
            "--vector-db-uri",
            envvar="MILVUS_URI",
            help="Milvus connection URI",
        ),
    ] = None,
    vector_db_token: Annotated[
        str | None,
        typer.Option(
            "--vector-db-token",
            envvar="MILVUS_TOKEN",
            help="Milvus authentication token",
        ),
    ] = None,
    vector_db_collection: Annotated[
        str,
        typer.Option(
            "--vector-db-collection",
            help="Collection name in VectorDB",
        ),
    ] = "pyeed",
    huggingface_token: Annotated[
        str | None,
        typer.Option(
            "--huggingface-token",
            envvar="HUGGINGFACE_TOKEN",
            help="HuggingFace API token (for gated models)",
        ),
    ] = None,
    model_name: Annotated[
        str,
        typer.Option(
            "--model-name",
            help="ESM2 model name from HuggingFace",
        ),
    ] = "facebook/esm2_t33_650M_UR50D",
    n_gpus: Annotated[
        int,
        typer.Option(
            "--n-gpus",
            help="Number of GPUs to use (-1 for all available)",
        ),
    ] = -1,
    model_dtype: Annotated[
        Literal["float16", "float32"],
        typer.Option(
            "--model-dtype",
            help="Model computation dtype",
        ),
    ] = "float32",
    return_dtype: Annotated[
        Literal["float16", "float32"],
        typer.Option(
            "--return-dtype",
            help="Embedding return dtype",
        ),
    ] = "float32",
    embedding_batch_size: Annotated[
        int,
        typer.Option(
            "--embedding-batch-size",
            help="GPU batch size for embedding",
        ),
    ] = 32,
    embedding_chunk_size: Annotated[
        int,
        typer.Option(
            "--embedding-chunk-size",
            help="Records to accumulate before embedding",
        ),
    ] = 1000,
    enrichment_batch_size: Annotated[
        int,
        typer.Option(
            "--enrichment-batch-size",
            help="Records to accumulate before enriching",
        ),
    ] = 1000,
    enrichment_max_concurrent: Annotated[
        int,
        typer.Option(
            "--enrichment-max-concurrent",
            help="Maximum concurrent API requests for enrichment",
        ),
    ] = 50,
) -> None:
    """Ingest proteins from UniProt by accession IDs.

    Fetches protein data from UniProt REST API, embeds sequences using ESM2,
    and enriches with taxonomy, reactions, and molecules.
    """
    ingest_uniprot(
        ids=ids,
        graph_db_uri=graph_db_uri,
        graph_db_user=graph_db_user,
        graph_db_password=graph_db_password,
        vector_db_uri=vector_db_uri,
        vector_db_token=vector_db_token,
        vector_db_collection=vector_db_collection,
        huggingface_token=huggingface_token,
        huggingface_model_name=model_name,
        huggingface_model_max_seq_length=1024,
        model_dtype=model_dtype,
        return_dtype=return_dtype,
        n_gpus=n_gpus,
        pooling_methods=None,
        embedding_chunk_size=embedding_chunk_size,
        embedding_batch_size=embedding_batch_size,
        uniprot_chunk_size=uniprot_chunk_size,
        uniprot_page_size=uniprot_page_size,
        enrichment_batch_size=enrichment_batch_size,
        enrichment_max_concurrent=enrichment_max_concurrent,
    )


@app.command(name="ingest-interpro")
def ingest_interpro_cmd(
    interpro_id: Annotated[
        str,
        typer.Argument(
            help="InterPro ID (e.g., IPR002133)",
        ),
    ],
    uniprot_chunk_size: Annotated[
        int,
        typer.Option(
            "--uniprot-chunk-size",
            help="Accessions per API batch",
        ),
    ] = 50,
    uniprot_page_size: Annotated[
        int,
        typer.Option(
            "--uniprot-page-size",
            help="Results per page",
        ),
    ] = 50,
    graph_db_uri: Annotated[
        str | None,
        typer.Option(
            "--graph-db-uri",
            envvar="NEO4J_URI",
            help="Neo4j connection URI",
        ),
    ] = None,
    graph_db_user: Annotated[
        str | None,
        typer.Option(
            "--graph-db-user",
            envvar="NEO4J_USER",
            help="Neo4j username",
        ),
    ] = None,
    graph_db_password: Annotated[
        str | None,
        typer.Option(
            "--graph-db-password",
            envvar="NEO4J_PASSWORD",
            help="Neo4j password",
        ),
    ] = None,
    vector_db_uri: Annotated[
        str | None,
        typer.Option(
            "--vector-db-uri",
            envvar="MILVUS_URI",
            help="Milvus connection URI",
        ),
    ] = None,
    vector_db_token: Annotated[
        str | None,
        typer.Option(
            "--vector-db-token",
            envvar="MILVUS_TOKEN",
            help="Milvus authentication token",
        ),
    ] = None,
    vector_db_collection: Annotated[
        str,
        typer.Option(
            "--vector-db-collection",
            help="Collection name in VectorDB",
        ),
    ] = "pyeed",
    huggingface_token: Annotated[
        str | None,
        typer.Option(
            "--huggingface-token",
            envvar="HUGGINGFACE_TOKEN",
            help="HuggingFace API token (for gated models)",
        ),
    ] = None,
    model_name: Annotated[
        str,
        typer.Option(
            "--model-name",
            help="ESM2 model name from HuggingFace",
        ),
    ] = "facebook/esm2_t33_650M_UR50D",
    n_gpus: Annotated[
        int,
        typer.Option(
            "--n-gpus",
            help="Number of GPUs to use (-1 for all available)",
        ),
    ] = -1,
    model_dtype: Annotated[
        Literal["float16", "float32"],
        typer.Option(
            "--model-dtype",
            help="Model computation dtype",
        ),
    ] = "float32",
    return_dtype: Annotated[
        Literal["float16", "float32"],
        typer.Option(
            "--return-dtype",
            help="Embedding return dtype",
        ),
    ] = "float32",
    embedding_batch_size: Annotated[
        int,
        typer.Option(
            "--embedding-batch-size",
            help="GPU batch size for embedding",
        ),
    ] = 32,
    embedding_chunk_size: Annotated[
        int,
        typer.Option(
            "--embedding-chunk-size",
            help="Records to accumulate before embedding",
        ),
    ] = 1000,
    enrichment_batch_size: Annotated[
        int,
        typer.Option(
            "--enrichment-batch-size",
            help="Records to accumulate before enriching",
        ),
    ] = 1000,
    enrichment_max_concurrent: Annotated[
        int,
        typer.Option(
            "--enrichment-max-concurrent",
            help="Maximum concurrent API requests for enrichment",
        ),
    ] = 50,
) -> None:
    """Ingest all proteins for an InterPro ID.

    Queries UniProt SPARQL endpoint to get all accessions for an InterPro family,
    then fetches full protein data, embeds sequences, and enriches with
    taxonomy, reactions, and molecules.
    """
    ingest_interpro(
        id=interpro_id,
        graph_db_uri=graph_db_uri,
        graph_db_user=graph_db_user,
        graph_db_password=graph_db_password,
        vector_db_uri=vector_db_uri,
        vector_db_token=vector_db_token,
        vector_db_collection=vector_db_collection,
        huggingface_token=huggingface_token,
        huggingface_model_name=model_name,
        huggingface_model_max_seq_length=1024,
        model_dtype=model_dtype,
        return_dtype=return_dtype,
        n_gpus=n_gpus,
        pooling_methods=None,
        embedding_chunk_size=embedding_chunk_size,
        embedding_batch_size=embedding_batch_size,
        uniprot_chunk_size=uniprot_chunk_size,
        uniprot_page_size=uniprot_page_size,
        enrichment_batch_size=enrichment_batch_size,
        enrichment_max_concurrent=enrichment_max_concurrent,
    )


if __name__ == "__main__":
    app()
