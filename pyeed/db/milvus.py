from __future__ import annotations

import os

import dotenv
import numpy as np
from loguru import logger
from pymilvus import (
    AsyncMilvusClient,
    DataType,
)

MILVUS_TO_NP_DTYPE_MAP = {
    DataType.FLOAT16_VECTOR: np.float16,
    DataType.FLOAT_VECTOR: np.float32,
}


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


async def initialize_collection_from_dict(
    collection_name: str,
    primary_field_name: str,
    record: dict[str, np.ndarray | str | int | float],
    *,
    client: AsyncMilvusClient,
) -> None:
    """Initialize a Milvus collection by inferring schema from a sample record.

    Creates a new collection with schema inferred from the provided record dictionary.
    The record serves as a template to determine field names, data types, and dimensions.
    For vector fields, automatically creates DISKANN indexes with COSINE similarity.

    Supported field types:
        - numpy.ndarray (float16/float32): Creates FLOAT16_VECTOR or FLOAT_VECTOR field
          with dimension inferred from array shape. Also creates a boolean `has_{field_name}`
          field to track presence of the vector.
        - str: Creates VARCHAR field with max_length=65535
        - int: Creates INT32 field (note: Python ints > INT32_MAX will overflow)
        - float: Creates FLOAT32 field

    Args:
        collection_name: Name of the collection to create.
        primary_field_name: Name of the primary key field (must exist in record).
        record: Sample record dictionary used to infer schema. Must contain
            primary_field_name and at least one other field. All fields in this record
            will be added to the schema.
        client: AsyncMilvusClient instance.

    Example:
        >>> client = AsyncMilvusClient(uri="...", token="...")
        >>> sample = {
        ...     "id": "protein_123",
        ...     "embedding": np.array([0.1, 0.2, 0.3], dtype=np.float32),
        ...     "name": "MyProtein",
        ...     "length": 100
        ... }
        >>> await initialize_collection_from_dict(
        ...     client, "proteins", "id", sample
        ... )
    """
    # Validate collection doesn't exist
    if collection_name in await client.list_collections():
        raise ValueError(f"Collection '{collection_name}' already exists")

    # Validate record is non-empty
    if not record:
        raise ValueError("Record cannot be empty")

    # Validate primary field exists
    if primary_field_name not in record:
        raise ValueError(
            f"Primary field '{primary_field_name}' not found in record. "
            f"Available fields: {list(record.keys())}"
        )

    schema = client.create_schema(
        auto_id=False,
        enable_dynamic_field=False,
    )

    index_params = client.prepare_index_params()

    # Add primary field (VARCHAR with max_length=128 for IDs)
    schema.add_field(
        field_name=primary_field_name,
        datatype=DataType.VARCHAR,
        max_length=128,
        is_primary=True,
    )

    # Infer field names, dtypes, and dimensions from record
    for key, value in record.items():
        if key == primary_field_name:
            continue

        if isinstance(value, np.ndarray):
            # Validate array is 1D
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

            schema.add_field(
                field_name=key,
                datatype=dtype,
                dim=dim,
            )
            index_params.add_index(
                field_name=key,
                index_name=f"{key}_index",
                index_type="DISKANN",
                metric_type="COSINE",
            )

        elif isinstance(value, str):
            dtype = DataType.VARCHAR
            max_length = 65535
            schema.add_field(
                field_name=key,
                datatype=dtype,
                max_length=max_length,
            )

        elif isinstance(value, int):
            # Note: Python int can exceed INT32 range, but Milvus INT32 is used
            # for compatibility. Consider INT64 if large values are expected.
            dtype = DataType.INT32
            schema.add_field(
                field_name=key,
                datatype=dtype,
            )

        elif isinstance(value, float):
            dtype = DataType.FLOAT32
            schema.add_field(
                field_name=key,
                datatype=dtype,
            )

        elif value is None:
            raise ValueError(
                f"Field '{key}' has None value. Cannot infer type from None. "
                f"Provide a sample value with the expected type."
            )

        else:
            raise ValueError(
                f"Unsupported dtype for field '{key}': {type(value)}. "
                f"Supported types: numpy.ndarray, str, int, float"
            )

    await client.create_collection(
        collection_name=collection_name,
        schema=schema,
        index_params=index_params,
        num_shards=8,
    )

    logger.info(f"Collection '{collection_name}' successfully created")


async def check_existing_ids(
    collection_name: str,
    ids: list[str],
    *,
    client: AsyncMilvusClient,
) -> set[str]:
    """Check if a collection contains the given ids.

    Args:
        collection_name: Name of the collection to check.
        ids: List of primary keys to check.
        client: AsyncMilvusClient instance.

    Returns:
        Set of primary keys that already exist in the collection.
    """
    schema = await client.describe_collection(collection_name)
    print(schema)
    response = await client.get(collection_name, ids=ids, output_fields=["protein_id"])
    return set(result["protein_id"] for result in response)


async def insert(
    collection_name: str,
    records: list[dict[str, int | float | str | np.ndarray]],
    *,
    client: AsyncMilvusClient,
) -> None:
    """Insert records into Milvus collection asynchronously.

    Automatically initializes collection if it doesn't exist.
    Converts records to Milvus format and inserts using async client.

    Args:
        collection_name: Name of the collection to insert into.
        records: List of records to insert.
        client: AsyncMilvusClient instance.
    """
    if not records:
        raise AssertionError("Records cannot be empty")

    await client.insert(collection_name=collection_name, data=records)
    logger.debug(f"Inserted {len(records)} records into '{collection_name}'")

    # # Search
    # # ------------------------------------------------------------

    # def get_vectors(
    #     self, collection_name: str, protein_ids: list[str], vector_field_name: str = "mean_pooling"
    # ) -> list[np.ndarray]:
    #     response = self.client.get(
    #         collection_name, ids=protein_ids, output_fields=[vector_field_name]
    #     )

    #     np_array = np.vstack([result[vector_field_name] for result in response]).astype(
    #         self._get_numpy_type(collection_name, vector_field_name)
    #     )
    #     return np_array

    # def vector_search(
    #     self,
    #     collection_name: str,
    #     query_vectors: list[np.ndarray],
    #     vector_field_name: str = "mean_pooling",
    #     return_vector: bool = False,
    #     n_hits: int = 10,
    #     return_fields: list[str] = [],
    #     radius: float = 1.0,
    #     range_filter: float = 0.0,
    #     hitlist_size: int = 16384,
    #     offset: int = 0,
    # ) -> list[list[dict[str, object]]]:
    #     if return_vector:
    #         return_fields.append(vector_field_name)

    #     return_fields = list(set(return_fields))

    #     # if vector is 1d make it 2d
    #     if len(query_vectors.shape) == 1:
    #         query_vectors = [query_vectors]

    #     records: list[dict]

    #     records = []
    #     records.extend(
    #         self.client.search(
    #             collection_name=collection_name,
    #             data=query_vectors,
    #             anns_field=vector_field_name,
    #             limit=n_hits,
    #             output_fields=return_fields,
    #             params={
    #                 "offset": offset,
    #                 "radius": radius,
    #                 "range_filter": range_filter,
    #             },
    #         )
    #     )

    #     # flatten the records
    #     clean_recs = []
    #     for rec in records[0]:
    #         value = rec["entity"][vector_field_name]
    #         vector = {
    #             vector_field_name: self._get_numpy_type(collection_name, vector_field_name)(value)
    #         }
    #         del rec["entity"]
    #         rec.update(vector)
    #         clean_recs.append(rec)

    #     df = pd.DataFrame.from_records(clean_recs)
    #     return df

    # def id_search(self, collection_name: str, protein_ids: list[str]) -> list[dict[str, object]]:
    #     if isinstance(protein_ids, str):
    #         protein_ids = [protein_ids]

    #     # get

    # # def _get_search_load_params(
    # #     self,
    # #     n_hits: int,
    # #     n_queries: int | None = None,
    # #     *,
    # #     max_total_hits_per_request: int = 16384,
    # #     request_size_limit_mb: float = 64.0,
    # #     safety: float = 0.6,
    # #     id_bytes: int = 8,
    # #     dist_bytes: int = 4,
    # #     wire_overhead: float = 1.6,
    # # ) -> tuple[int, int]:
    # #     """Return (hitlist_size, queries_per_search) without exceeding limits.

    # #     n_hits: desired results per query (topK).
    # #     n_queries: optional cap on queries per network request.
    # #     """

    # #     # Byte-budget-derived cap on total results we can return in one request.
    # #     request_budget_bytes = int(request_size_limit_mb * (1024**2) * safety)
    # #     bytes_per_result = int((id_bytes + dist_bytes) * wire_overhead)
    # #     budget_cap_total_hits = request_budget_bytes // bytes_per_result

    # #     # Effective total-results cap per request respects both byte budget and server max.
    # #     total_hits_cap = min(budget_cap_total_hits, max_total_hits_per_request)

    # #     # Per-query result count cannot exceed requested n_hits nor the total cap.
    # #     hitlist_size = min(n_hits, total_hits_cap)

    # #     # Number of queries we can pack while staying within the total results cap.
    # #     queries_per_search = total_hits_cap // hitlist_size

    # #     if n_queries is not None:
    # #         queries_per_search = min(queries_per_search, n_queries)

    # #     # Ensure product stays within the cap (guard against rounding).
    # #     if queries_per_search * hitlist_size > total_hits_cap:
    # #         queries_per_search = total_hits_cap // hitlist_size

    # #     return hitlist_size, queries_per_search

    # def get_all_vectors(
    #     self,
    #     collection_name: str,
    #     vector_field_name: str = "mean_pooling",
    #     show_progress: bool = True,
    # ) -> list[np.ndarray]:
    #     it = self.client.query_iterator(
    #         collection_name=collection_name,
    #         batch_size=1000,
    #         output_fields=[vector_field_name],
    #     )
    #     protein_ids = []
    #     vectors = []

    #     # get total number of rows
    #     total_rows = self.client.get_collection_stats(collection_name)["row_count"]

    #     with create_progress() as progress:
    #         task = progress.add_task(
    #             f"Loading {vector_field_name} vectors from {collection_name}", total=total_rows
    #         )

    #         while True:
    #             batch = it.next()
    #             if not batch:
    #                 break

    #             # extract protein id and make in list of protein ids
    #             for rec in batch:
    #                 protein_ids.append(rec["protein_id"])
    #                 vectors.append(rec[vector_field_name])

    #             progress.update(task, advance=len(batch))

    #         task_convert = progress.add_task(
    #             "Converting vectors to numpy array", total=None, start=True
    #         )
    #         vectors = np.vstack(vectors).astype(
    #             self._get_numpy_type(collection_name, vector_field_name)
    #         )
    #         progress.update(task_convert, completed=True)
    #         sleep(1)

    #         return protein_ids, vectors

    # def _get_numpy_type(
    #     self,
    #     collection_name: str,
    #     vector_field_name: str,
    # ) -> np.dtype:
    #     schema = self.client.describe_collection(collection_name)
    #     for field in schema["fields"]:
    #         if field["name"] == vector_field_name:
    #             return MILVUS_TO_NP_DTYPE_MAP[field["type"]]
    #     raise ValueError(
    #         f"Vector field '{vector_field_name}' not found in collection '{collection_name}'"
    #     )


if __name__ == "__main__":
    import asyncio

    from rich import print

    test_records = [
        {"id": "1", "embedding": np.array([0.1, 0.2, 0.3], dtype=np.float32)},
        {"id": "2", "embedding": np.array([0.4, 0.5, 0.6], dtype=np.float32)},
        {"id": "3", "embedding": np.array([0, 0, 0], dtype=np.float32)},
    ]

    collection_name = "test_v2"

    # test all functions
    async def main() -> None:
        print("Initializing collection")
        client = get_async_milvus_client()

        # drop collection if it exists
        await client.drop_collection(collection_name)
        print("Collection dropped")

        await initialize_collection_from_dict(
            collection_name=collection_name,
            primary_field_name="id",
            record=test_records[0],
            client=client,
        )
        print("Inserting test records")
        await insert(
            collection_name=collection_name,
            records=test_records,
            client=client,
        )

    asyncio.run(main())
