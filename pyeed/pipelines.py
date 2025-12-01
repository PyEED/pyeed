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
from pyeed.db.neo4j import get_async_driver
from pyeed.embed.esm2 import ESM2Embedder
from pyeed.embed.pooling import PoolingLike, mean_pooling

from .ingest.core.pipeline import Pipeline
from .ingest.sources.fasta import build_header_index
from .ingest.stages.embed import EmbeddingStage
from .ingest.stages.enrich_molecules import MoleculeEnrichmentStage
from .ingest.stages.enrich_reactions import ReactionEnrichmentStage
from .ingest.stages.enrich_taxonomy import TaxonomyEnrichmentStage
from .ingest.stages.fasta import FASTAReaderStage
from .ingest.stages.sink import MilvusUpsertStage, Neo4jUpsertStage
from .ingest.stages.uniprot import InterProReaderStage, UniProtReaderStage
from .utils.progress import CONSOLE


def ingest_fasta(
    fasta_path: str,
    chunk_size: int = 512,
    header_fn: Callable[[str], str] | None = None,
    taxon_fn: Callable[[str], str] | None = None,
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
    pooling_methods: PoolingLike | None = None,
    embedding_chunk_size: int = 1000,
    embedding_batch_size: int = 32,
    enrichment_batch_size: int = 50,
    enrichment_max_concurrent: int = 50,
) -> None:
    """Ingest a FASTA file into the pipeline.

    Args:
        fasta_path: Path to FASTA file
        chunk_size: Number of sequences per chunk for FASTA reading
        header_fn: Function to extract protein_id from header
            Signature: (header: str) -> str
            If None, uses entire header as id of protein
        taxon_fn: Function to extract taxon_id from header
            Signature: (header: str) -> str
            If None, taxon_id will be None
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
        enrichment_batch_size: Records to accumulate before enriching (default: 1000)
        enrichment_max_concurrent: Maximum concurrent API requests for enrichment (default: 10)

    Returns:
        None
    """
    if pooling_methods is None:
        pooling_methods = [mean_pooling]

    bar = Progress(
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        TimeElapsedColumn(),
        transient=False,
    )
    spin = Progress(
        SpinnerColumn(style="cyan", spinner_name="dots", finished_text="✓"),
        TextColumn("{task.description}"),
        transient=False,
    )

    # Tasks
    graph_db_connect_task = spin.add_task("Connecting to GraphDB", total=1)
    vector_db_connect_task = spin.add_task("Connecting to VectorDB", total=1)
    embedder_init_task = spin.add_task("Initializing ESM2Embedder", total=1)
    read_task = bar.add_task("Read FASTA", total=None)
    write_task = bar.add_task("Add Sequences to GraphDB", total=None)
    embed_task = bar.add_task("Embed Sequences", total=None)
    milvus_task = bar.add_task("Add Embeddings to VectorDB", total=None)
    enrich_taxonomy_task = bar.add_task("Enrich Taxonomy", total=None)

    with Live(Group(spin, bar), console=CONSOLE, refresh_per_second=8):
        # Connect to GraphDB
        driver = get_async_driver(
            uri=graph_db_uri,
            user=graph_db_user,
            password=graph_db_password,
        )
        logger.info(f"Connected to GraphDB: {graph_db_uri}")
        spin.update(graph_db_connect_task, completed=1)
        spin.update(graph_db_connect_task, description=f"Connected to GraphDB: {graph_db_uri}")

        # Connect to VectorDB
        vector_db = VectorDB(
            uri=vector_db_uri,
            token=vector_db_token,
            collection_name=vector_db_collection,
        )
        logger.info(f"Connected to VectorDB: {vector_db.uri}")
        spin.update(vector_db_connect_task, completed=1)
        spin.update(vector_db_connect_task, description=f"Connected to VectorDB: {vector_db.uri}")

        # Initialize ESM2Embedder (only if embedding is enabled)
        embedder = ESM2Embedder(
            model_name=huggingface_model_name,
            model_dtype=model_dtype,
            return_dtype=return_dtype,
            n_gpus=n_gpus,
            pooling_methods=pooling_methods,
            max_length=huggingface_model_max_seq_length,
            huggingface_token=huggingface_token,
        )
        asyncio.run(embedder.initialize())
        logger.info(f"Initialized {embedder.model_name} on {len(embedder.devices)} devices")
        spin.update(embedder_init_task, completed=1)
        spin.update(
            embedder_init_task,
            description=f"Initialized {embedder.model_name} on {len(embedder.devices)} devices",
        )

        # Build header index from FASTA file
        offsets = build_header_index(fasta_path)
        total_sequences = len(offsets)
        bar.update(read_task, total=total_sequences)
        bar.update(embed_task, total=total_sequences)
        bar.update(write_task, total=total_sequences)
        bar.update(milvus_task, total=total_sequences)
        bar.update(enrich_taxonomy_task, total=total_sequences)

        # Build pipeline
        pipeline = Pipeline(progress=bar)
        neo4j_queue = pipeline.add_queue("neo4j", maxsize=3000)
        embedding_queue = pipeline.add_queue("embedding", maxsize=3000)
        taxonomy_queue = pipeline.add_queue("taxonomy", maxsize=3000)
        milvus_queue = pipeline.add_queue("milvus", maxsize=3000)

        pipeline.add_stage(
            stage=FASTAReaderStage(
                fasta_path=fasta_path,
                offsets=offsets,
                chunk_size=chunk_size,
                header_extractor=header_fn,
                taxon_extractor=taxon_fn,
            ),
            input_queues=[],
            output_queues=[neo4j_queue],
            task_id=read_task,
        )

        pipeline.add_stage(
            stage=Neo4jUpsertStage(driver=driver, batch_size=512),
            input_queues=[neo4j_queue],
            output_queues=[embedding_queue, taxonomy_queue],
            task_id=write_task,
        )

        pipeline.add_stage(
            stage=EmbeddingStage(
                embedder=embedder,
                chunk_size=embedding_chunk_size,
                batch_size=embedding_batch_size,
                vector_db=vector_db,
                collection_name=vector_db_collection,
            ),
            input_queues=[embedding_queue],
            output_queues=[milvus_queue],
            task_id=embed_task,
        )

        pipeline.add_stage(
            stage=MilvusUpsertStage(
                vector_db=vector_db,
                collection_name=vector_db_collection,
                batch_size=512,
            ),
            input_queues=[milvus_queue],
            output_queues=[],
            task_id=milvus_task,
        )

        # Create shared DB semaphore for enrichment stages
        db_semaphore = asyncio.Semaphore(1)

        pipeline.add_stage(
            stage=TaxonomyEnrichmentStage(
                driver=driver,
                db_semaphore=db_semaphore,
                batch_size=enrichment_batch_size,
                max_concurrent=enrichment_max_concurrent,
            ),
            input_queues=[taxonomy_queue],
            output_queues=[],
            task_id=enrich_taxonomy_task,
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
    huggingface_token: str | None = None,
    huggingface_model_name: str = "facebook/esm2_t33_650M_UR50D",
    huggingface_model_max_seq_length: int = 1024,
    model_dtype: Literal["float16", "float32"] = "float32",
    return_dtype: Literal["float16", "float32"] = "float32",
    n_gpus: int = -1,
    pooling_methods: PoolingLike | None = None,
    embedding_chunk_size: int = 1000,
    embedding_batch_size: int = 32,
    uniprot_chunk_size: int = 50,
    uniprot_page_size: int = 50,
    enrichment_batch_size: int = 1000,
    enrichment_max_concurrent: int = 50,
) -> None:
    """Ingest proteins from UniProt by accession IDs.

    Args:
        ids: List of UniProt accession IDs (e.g., ["P12345", "Q9Y6X9"])
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
        uniprot_chunk_size: Accessions per API batch (default: 50)
        uniprot_page_size: Results per page (default: 50)
        enrichment_batch_size: Records to accumulate before enriching (default: 1000)
        enrichment_max_concurrent: Maximum concurrent API requests for enrichment (default: 10)

    Returns:
        None
    """
    if pooling_methods is None:
        pooling_methods = [mean_pooling]

    bar = Progress(
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        TimeElapsedColumn(),
        transient=False,
    )
    spin = Progress(
        SpinnerColumn(style="cyan", spinner_name="dots", finished_text="✓"),
        TextColumn("{task.description}"),
        transient=False,
    )

    # Tasks
    graph_db_connect_task = spin.add_task("Connecting to GraphDB", total=1)
    vector_db_connect_task = spin.add_task("Connecting to VectorDB", total=1)
    embedder_init_task = spin.add_task("Initializing ESM2Embedder", total=1)

    fetch_task = bar.add_task("Fetch from UniProt", total=len(ids))
    write_task = bar.add_task("Add Sequences to GraphDB", total=len(ids))
    embed_task = bar.add_task("Embed Sequences", total=len(ids))
    milvus_task = bar.add_task("Add Embeddings to VectorDB", total=len(ids))
    enrich_taxonomy_task = bar.add_task("Enrich Taxonomy", total=len(ids))
    enrich_reactions_task = bar.add_task("Enrich Reactions", total=None)
    enrich_molecules_task = bar.add_task("Enrich Molecules", total=None)

    with Live(Group(spin, bar), console=CONSOLE, refresh_per_second=8):
        # Connect to GraphDB
        driver = get_async_driver(
            uri=graph_db_uri,
            user=graph_db_user,
            password=graph_db_password,
        )
        logger.info(f"Connected to GraphDB: {graph_db_uri}")
        spin.update(graph_db_connect_task, completed=1)
        spin.update(graph_db_connect_task, description=f"Connected to GraphDB: {graph_db_uri}")

        # Connect to VectorDB
        vector_db = VectorDB(
            uri=vector_db_uri,
            token=vector_db_token,
            collection_name=vector_db_collection,
        )
        spin.update(vector_db_connect_task, completed=1)
        spin.update(vector_db_connect_task, description=f"Connected to VectorDB: {vector_db.uri}")
        logger.info(f"Connected to VectorDB: {vector_db.uri}")

        embedder = ESM2Embedder(
            model_name=huggingface_model_name,
            model_dtype=model_dtype,
            return_dtype=return_dtype,
            n_gpus=n_gpus,
            pooling_methods=pooling_methods,
            max_length=huggingface_model_max_seq_length,
            huggingface_token=huggingface_token,
        )
        asyncio.run(embedder.initialize())
        logger.info(f"Initialized {embedder.model_name} on {len(embedder.devices)} devices")
        spin.update(embedder_init_task, completed=1)
        spin.update(
            embedder_init_task,
            description=f"Initialized {embedder.model_name} on {len(embedder.devices)} devices",
        )

        # Build pipeline
        pipeline = Pipeline(progress=bar)

        neo4j_queue = pipeline.add_queue("neo4j", maxsize=1000)
        embedding_queue = pipeline.add_queue("embedding", maxsize=1000)
        taxonomy_queue = pipeline.add_queue("taxonomy", maxsize=1000)
        reaction_queue = pipeline.add_queue("reaction", maxsize=1000)
        molecule_queue = pipeline.add_queue("molecule", maxsize=1000)
        milvus_queue = pipeline.add_queue("milvus", maxsize=1000)

        # Stage 1: Fetch from UniProt
        pipeline.add_stage(
            stage=UniProtReaderStage(
                accessions=ids,
                chunk_size=uniprot_chunk_size,
                size_per_page=uniprot_page_size,
            ),
            input_queues=[],
            output_queues=[neo4j_queue],
            task_id=fetch_task,
        )

        # Stage 2: Neo4j (forwards to embedding, taxonomy, and reactions)
        pipeline.add_stage(
            stage=Neo4jUpsertStage(driver=driver, batch_size=100),
            input_queues=[neo4j_queue],
            output_queues=[embedding_queue, taxonomy_queue, reaction_queue],
            task_id=write_task,
        )

        # Path 1: Embedding pipeline
        # Stage 3a: Embedding (checks Milvus, embeds new)
        pipeline.add_stage(
            stage=EmbeddingStage(
                embedder=embedder,
                chunk_size=embedding_chunk_size,
                batch_size=embedding_batch_size,
                vector_db=vector_db,
                collection_name=vector_db_collection,
            ),
            input_queues=[embedding_queue],
            output_queues=[milvus_queue],
            task_id=embed_task,
        )

        # Stage 4a: Milvus (inserts embeddings)
        pipeline.add_stage(
            stage=MilvusUpsertStage(
                vector_db=vector_db,
                collection_name=vector_db_collection,
                batch_size=2000,
            ),
            input_queues=[milvus_queue],
            output_queues=[],
            task_id=milvus_task,
        )

        # Path 2: Enrichment pipeline (taxonomy → reactions → molecules)
        # Create shared DB semaphore for all enrichment stages
        db_semaphore = asyncio.Semaphore(1)

        pipeline.add_stage(
            stage=TaxonomyEnrichmentStage(
                driver=driver,
                db_semaphore=db_semaphore,
                batch_size=enrichment_batch_size,
                max_concurrent=enrichment_max_concurrent,
            ),
            input_queues=[taxonomy_queue],
            output_queues=[],
            task_id=enrich_taxonomy_task,
        )

        pipeline.add_stage(
            stage=ReactionEnrichmentStage(
                driver=driver,
                db_semaphore=db_semaphore,
                batch_size=enrichment_batch_size,
                max_concurrent=enrichment_max_concurrent,
            ),
            input_queues=[reaction_queue],
            output_queues=[molecule_queue],
            task_id=enrich_reactions_task,
        )

        pipeline.add_stage(
            stage=MoleculeEnrichmentStage(
                driver=driver,
                db_semaphore=db_semaphore,
                batch_size=enrichment_batch_size,
                max_concurrent=enrichment_max_concurrent,
            ),
            input_queues=[molecule_queue],
            output_queues=[],
            task_id=enrich_molecules_task,
        )

        asyncio.run(pipeline.run())


def ingest_interpro(
    id: str,
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
    pooling_methods: PoolingLike | None = None,
    embedding_chunk_size: int = 1000,
    embedding_batch_size: int = 32,
    uniprot_chunk_size: int = 50,
    uniprot_page_size: int = 50,
    enrichment_batch_size: int = 1000,
    enrichment_max_concurrent: int = 50,
) -> None:
    """Ingest all proteins for an InterPro ID.

    First queries the UniProt SPARQL endpoint to get all accession IDs
    associated with an InterPro family, then fetches full protein data
    for each accession.

    Args:
        id: InterPro ID (e.g., "IPR002133")
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
        uniprot_chunk_size: Accessions per API batch (default: 50)
        uniprot_page_size: Results per page (default: 50)
        enrichment_batch_size: Records to accumulate before enriching (default: 1000)
        enrichment_max_concurrent: Maximum concurrent API requests for enrichment (default: 10)

    Returns:
        None
    """
    if pooling_methods is None:
        pooling_methods = [mean_pooling]

    bar = Progress(
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        TimeElapsedColumn(),
        transient=False,
    )
    spin = Progress(
        SpinnerColumn(style="cyan", spinner_name="dots", finished_text="✓"),
        TextColumn("{task.description}"),
        transient=False,
    )

    # Tasks
    graph_db_connect_task = spin.add_task("Connecting to GraphDB", total=1)
    vector_db_connect_task = spin.add_task("Connecting to VectorDB", total=1)
    embedder_init_task = spin.add_task("Initializing ESM2Embedder", total=1)

    fetch_task = bar.add_task(f"Fetch from InterPro ({id})", total=None)
    write_task = bar.add_task("Add Sequences to GraphDB", total=None)
    embed_task = bar.add_task("Embed Sequences", total=None)
    milvus_task = bar.add_task("Add Embeddings to VectorDB", total=None)
    enrich_taxonomy_task = bar.add_task("Enrich Taxonomy", total=None)
    enrich_reactions_task = bar.add_task("Enrich Reactions", total=None)
    enrich_molecules_task = bar.add_task("Enrich Molecules", total=None)

    with Live(Group(spin, bar), console=CONSOLE, refresh_per_second=8):
        # Connect to GraphDB
        driver = get_async_driver(
            uri=graph_db_uri,
            user=graph_db_user,
            password=graph_db_password,
        )
        logger.info(f"Connected to GraphDB: {graph_db_uri}")
        spin.update(graph_db_connect_task, completed=1)
        spin.update(graph_db_connect_task, description=f"Connected to GraphDB: {graph_db_uri}")

        # Connect to VectorDB
        vector_db = VectorDB(
            uri=vector_db_uri,
            token=vector_db_token,
            collection_name=vector_db_collection,
        )
        logger.info(f"Connected to VectorDB: {vector_db.uri}")
        spin.update(vector_db_connect_task, completed=1)
        spin.update(vector_db_connect_task, description=f"Connected to VectorDB: {vector_db.uri}")

        # Initialize ESM2Embedder (only if embedding is enabled)

        embedder = ESM2Embedder(
            model_name=huggingface_model_name,
            model_dtype=model_dtype,
            return_dtype=return_dtype,
            n_gpus=n_gpus,
            pooling_methods=pooling_methods,
            max_length=huggingface_model_max_seq_length,
            huggingface_token=huggingface_token,
        )
        asyncio.run(embedder.initialize())
        logger.info(f"Initialized {embedder.model_name} on {len(embedder.devices)} devices")
        spin.update(embedder_init_task, completed=1)
        spin.update(
            embedder_init_task,
            description=f"Initialized {embedder.model_name} on {len(embedder.devices)} devices",
        )

        # Build pipeline
        pipeline = Pipeline(progress=bar)

        neo4j_queue = pipeline.add_queue("neo4j", maxsize=1400)
        embedding_queue = pipeline.add_queue("embedding", maxsize=2000)
        taxonomy_queue = pipeline.add_queue("taxonomy", maxsize=1500)
        reaction_queue = pipeline.add_queue("reaction", maxsize=2000)
        molecule_queue = pipeline.add_queue("molecule", maxsize=3000)
        milvus_queue = pipeline.add_queue("milvus", maxsize=4000)

        pipeline.add_stage(
            stage=InterProReaderStage(
                interpro_id=id,
                chunk_size=uniprot_chunk_size,
                size_per_page=uniprot_page_size,
            ),
            input_queues=[],
            output_queues=[neo4j_queue],
            task_id=fetch_task,
        )

        pipeline.add_stage(
            stage=Neo4jUpsertStage(driver=driver, batch_size=100),
            input_queues=[neo4j_queue],
            output_queues=[embedding_queue, taxonomy_queue, reaction_queue],
            task_id=write_task,
        )

        pipeline.add_stage(
            stage=EmbeddingStage(
                embedder=embedder,
                chunk_size=embedding_chunk_size,
                batch_size=embedding_batch_size,
                vector_db=vector_db,
                collection_name=vector_db_collection,
            ),
            input_queues=[embedding_queue],
            output_queues=[milvus_queue],
            task_id=embed_task,
        )

        pipeline.add_stage(
            stage=MilvusUpsertStage(
                vector_db=vector_db,
                collection_name=vector_db_collection,
                batch_size=2000,
            ),
            input_queues=[milvus_queue],
            output_queues=[],
            task_id=milvus_task,
        )

        # Create shared DB semaphore for all enrichment stages
        db_semaphore = asyncio.Semaphore(1)

        pipeline.add_stage(
            stage=TaxonomyEnrichmentStage(
                driver=driver,
                db_semaphore=db_semaphore,
                batch_size=enrichment_batch_size,
                max_concurrent=enrichment_max_concurrent,
            ),
            input_queues=[taxonomy_queue],
            output_queues=[],
            task_id=enrich_taxonomy_task,
        )

        pipeline.add_stage(
            stage=ReactionEnrichmentStage(
                driver=driver,
                db_semaphore=db_semaphore,
                batch_size=enrichment_batch_size,
                max_concurrent=enrichment_max_concurrent,
            ),
            input_queues=[reaction_queue],
            output_queues=[molecule_queue],
            task_id=enrich_reactions_task,
        )

        pipeline.add_stage(
            stage=MoleculeEnrichmentStage(
                driver=driver,
                db_semaphore=db_semaphore,
                batch_size=enrichment_batch_size,
                max_concurrent=enrichment_max_concurrent,
            ),
            input_queues=[molecule_queue],
            output_queues=[],
            task_id=enrich_molecules_task,
        )

        asyncio.run(pipeline.run())


if __name__ == "__main__":

    def extract_rep_id(header: str) -> str:
        """
        From:
            'UniRef90_UPI002E2621C6 uncharacterized protein ...'
        → 'UPI002E2621C6'
        """
        return header.split("_", 1)[1].split(" ", 1)[0]

    def extract_taxid(header: str) -> str:
        """
        From:
            '... Tax=Corticium candelabrum TaxID=121492 RepID=...'
        → '121492'
        """
        return header.split("TaxID=", 1)[1].split()[0]

    fasta_path = "/home/mha/downloads/uniref_parts/uniref90_part03.fasta"
    ingest_fasta(
        fasta_path=fasta_path,
        header_fn=extract_rep_id,
        taxon_fn=extract_taxid,
        n_gpus=2,
    )
