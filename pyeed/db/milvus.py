from __future__ import annotations

import os
from typing import Literal

import dotenv
import numpy as np
from loguru import logger
from pymilvus import AsyncMilvusClient, DataType

MILVUS_TO_NP_DTYPE_MAP = {
    DataType.FLOAT16_VECTOR: np.float16,
    DataType.FLOAT_VECTOR: np.float32,
}

type IndexType = Literal["DISKANN", "HNSW", "IVF_FLAT", "FLAT"]
type MetricType = Literal["COSINE", "L2", "IP"]


def get_async_milvus_client(
    uri: str | None = None,
    token: str | None = None,
) -> AsyncMilvusClient:
    """Create and return an AsyncMilvusClient instance.

    Loads connection parameters from environment variables if not provided.
    Environment variables are loaded from .env file if present.

    Args:
        uri: Milvus server URI. If None, reads from MILVUS_URL environment variable.
        token: Milvus authentication token. If None, reads from MILVUS_TOKEN
            environment variable.

    Returns:
        AsyncMilvusClient instance configured with the provided or environment
        URI and token.
    """
    dotenv.load_dotenv()
    uri = uri or os.getenv("MILVUS_URL")
    token = token or os.getenv("MILVUS_TOKEN")
    client = AsyncMilvusClient(uri=uri, token=token)
    logger.info("Connected to Milvus", extra={"uri": uri})
    return client


async def create_collection_for_bulk_load(
    collection_name: str,
    primary_field_name: str,
    record: dict[str, np.ndarray | str | int | float],
    *,
    client: AsyncMilvusClient,
) -> None:
    """Create a Milvus collection WITHOUT indexes for bulk loading.

    Creates a new collection with schema inferred from the provided record dictionary.
    No indexes are created - call `create_vector_index()` after bulk loading is complete.

    Supported field types:
        - numpy.ndarray (float16/float32): Creates FLOAT16_VECTOR or FLOAT_VECTOR field
        - str: Creates VARCHAR field with max_length=65535
        - int: Creates INT32 field
        - float: Creates FLOAT32 field

    Args:
        collection_name: Name of the collection to create.
        primary_field_name: Name of the primary key field (must exist in record).
        record: Sample record dictionary used to infer schema.
        client: AsyncMilvusClient instance.

    Example:
        >>> sample = {"id": "protein_123", "embedding": np.zeros(1280, dtype=np.float16)}
        >>> await create_collection_for_bulk_load(
        ...     "proteins", "id", sample, client=client
        ... )
        >>> # ... bulk insert data ...
        >>> await create_vector_index("proteins", "embedding", client=client)
    """
    if collection_name in await client.list_collections():
        raise ValueError(f"Collection '{collection_name}' already exists")

    if not record:
        raise ValueError("Record cannot be empty")

    if primary_field_name not in record:
        raise ValueError(
            f"Primary field '{primary_field_name}' not found in record. "
            f"Available fields: {list(record.keys())}"
        )

    schema = client.create_schema(auto_id=False, enable_dynamic_field=False)

    # Add primary field
    schema.add_field(
        field_name=primary_field_name,
        datatype=DataType.VARCHAR,
        max_length=128,
        is_primary=True,
    )

    # Infer field types from record
    for key, value in record.items():
        if key == primary_field_name:
            continue

        if isinstance(value, np.ndarray):
            if value.ndim != 1:
                raise ValueError(f"Vector field '{key}' must be 1D array, got shape {value.shape}")

            dim = value.shape[0]
            if value.dtype == np.float16:
                dtype = DataType.FLOAT16_VECTOR
            elif value.dtype == np.float32:
                dtype = DataType.FLOAT_VECTOR
            else:
                raise ValueError(
                    f"Unsupported vector dtype for field '{key}': {value.dtype}. "
                    f"Expected float16 or float32"
                )
            schema.add_field(field_name=key, datatype=dtype, dim=dim)

        elif isinstance(value, str):
            schema.add_field(field_name=key, datatype=DataType.VARCHAR, max_length=65535)

        elif isinstance(value, int):
            schema.add_field(field_name=key, datatype=DataType.INT32)

        elif isinstance(value, float):
            schema.add_field(field_name=key, datatype=DataType.FLOAT32)

        elif value is None:
            raise ValueError(f"Field '{key}' has None value. Cannot infer type from None.")

        else:
            raise ValueError(
                f"Unsupported dtype for field '{key}': {type(value)}. "
                f"Supported types: numpy.ndarray, str, int, float"
            )

    # Create collection WITHOUT index for bulk loading
    await client.create_collection(
        collection_name=collection_name,
        schema=schema,
        num_shards=8,
    )

    logger.info(f"Collection '{collection_name}' created (no index, ready for bulk load)")


async def create_vector_index(
    collection_name: str,
    vector_field_name: str,
    *,
    client: AsyncMilvusClient,
    index_type: IndexType = "DISKANN",
    metric_type: MetricType = "COSINE",
) -> None:
    """Create a vector index on a collection after bulk loading.

    Call this after all data has been inserted to build the index once.
    Building the index after bulk load is much faster than incremental indexing.

    Args:
        collection_name: Name of the collection.
        vector_field_name: Name of the vector field to index.
        client: AsyncMilvusClient instance.
        index_type: Type of index (DISKANN, HNSW, IVF_FLAT, FLAT).
        metric_type: Distance metric (COSINE, L2, IP).

    Example:
        >>> await create_vector_index("proteins", "embedding", client=client)
    """
    index_params = client.prepare_index_params()
    index_params.add_index(
        field_name=vector_field_name,
        index_name=f"{vector_field_name}_index",
        index_type=index_type,
        metric_type=metric_type,
    )

    await client.create_index(collection_name=collection_name, index_params=index_params)
    logger.info(
        f"Created {index_type} index on '{collection_name}.{vector_field_name}'",
        extra={"metric": metric_type},
    )


async def insert_batch(
    collection_name: str,
    ids: list[str],
    embeddings: np.ndarray,
    *,
    client: AsyncMilvusClient,
    id_field: str = "id",
    embedding_field: str = "embedding",
) -> None:
    """Insert batch using column-oriented format for maximum throughput.

    Args:
        collection_name: Name of the collection to insert into.
        ids: List of primary key values.
        embeddings: 2D numpy array of shape (n, dim), dtype float16 or float32.
            Must be contiguous in memory (C-order).
        client: AsyncMilvusClient instance.
        id_field: Name of the primary key field.
        embedding_field: Name of the embedding field.

    Example:
        >>> ids = ["p1", "p2", "p3"]
        >>> embeddings = np.random.randn(3, 1280).astype(np.float16)
        >>> await insert_batch("proteins", ids, embeddings, client=client)
    """
    if len(ids) != embeddings.shape[0]:
        raise ValueError(f"Mismatch: {len(ids)} ids vs {embeddings.shape[0]} embeddings")

    # Ensure contiguous memory layout
    if not embeddings.flags["C_CONTIGUOUS"]:
        embeddings = np.ascontiguousarray(embeddings)

    import time

    start_time = time.perf_counter()
    # Make data a list of dicts (row-oriented), each with id and embedding
    data = [
        {id_field: id_, embedding_field: emb} for id_, emb in zip(ids, embeddings, strict=False)
    ]
    elapsed_time = time.perf_counter() - start_time
    print(f"Batch data preparation to dict took {elapsed_time:.6f} seconds")

    await client.insert(collection_name=collection_name, data=data)
    logger.debug(f"Inserted {len(ids)} records into '{collection_name}'")


async def check_existing_ids(
    collection_name: str,
    ids: list[str],
    *,
    client: AsyncMilvusClient,
    id_field: str = "id",
) -> set[str]:
    """Check if a collection contains the given ids.

    Args:
        collection_name: Name of the collection to check.
        ids: List of primary keys to check.
        client: AsyncMilvusClient instance.
        id_field: Name of the primary key field.

    Returns:
        Set of primary keys that already exist in the collection.
    """
    response = await client.get(collection_name, ids=ids, output_fields=[id_field])
    return set(result[id_field] for result in response)


async def insert(
    collection_name: str,
    records: list[dict[str, int | float | str | np.ndarray]],
    *,
    client: AsyncMilvusClient,
) -> None:
    """Insert records into Milvus collection (row-oriented format).

    For bulk operations, prefer `bulk_insert()` with column-oriented data.

    Args:
        collection_name: Name of the collection to insert into.
        records: List of records to insert.
        client: AsyncMilvusClient instance.
    """
    if not records:
        raise AssertionError("Records cannot be empty")

    await client.insert(collection_name=collection_name, data=records)
    logger.debug(f"Inserted {len(records)} records into '{collection_name}'")


# Legacy alias for backward compatibility
initialize_collection_from_dict = create_collection_for_bulk_load
