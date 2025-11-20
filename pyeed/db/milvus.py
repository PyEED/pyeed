from __future__ import annotations

import asyncio
import os
from dataclasses import dataclass
from time import sleep

import dotenv
import numpy as np
import pandas as pd
from loguru import logger
from pymilvus import (
    AsyncMilvusClient,
    Collection,
    DataType,
    MilvusClient,
)

from pyeed.utils.progress import create_progress

from ..embed.types import EmbeddingBatch

MILVUS_TO_NP_DTYPE_MAP = {
    DataType.FLOAT16_VECTOR: np.float16,
    DataType.FLOAT_VECTOR: np.float32,
}


@dataclass
class UploadStats:
    """Statistics for upload operations."""

    total_inserted: int = 0
    total_batches: int = 0
    failed_inserts: int = 0


class VectorDB:
    MAX_HITLIST_SIZE = 16384

    def __init__(
        self,
        uri: str | None = None,
        token: str | None = None,
        collection_name: str = "pyeed",
    ):
        dotenv.load_dotenv()
        uri = uri or os.getenv("MILVUS_URI")
        token = token or os.getenv("MILVUS_TOKEN")

        self.uri = uri
        self.collection_name = collection_name

        self.async_client, self.client = self._connect(uri, token)
        self._connected = True

        self.collections = self.client.list_collections()
        self.databases = self.client.list_databases()
        self._initialized_collections: set[str] = set()

    def _connect(self, uri: str | None, token: str | None):
        if not uri or not token:
            raise ValueError("Both 'uri' and 'token' must be provided to connect to Milvus.")

        async_client = AsyncMilvusClient(uri=uri, token=token)
        client = MilvusClient(uri=uri, token=token)
        logger.info("Connected to Milvus", extra={"uri": uri, "token": token})
        return async_client, client

    def create_collection(
        self,
        collection_name: str,
        vec_field_names: list[str],
        vec_field_dtypes: list[DataType],
        vec_dims: list[int],
        include_sequence: bool,
    ) -> Collection:
        schema = self.client.create_schema(
            auto_id=False,
            enable_dynamic_field=False,
        )
        index_params = self.client.prepare_index_params()

        schema.add_field(
            field_name="protein_id",
            datatype=DataType.VARCHAR,
            max_length=65535,
            is_primary=True,
        )

        for name, dtype, dim in zip(vec_field_names, vec_field_dtypes, vec_dims, strict=True):
            schema.add_field(
                field_name=name,
                datatype=dtype,
                dim=dim,
            )
            index_params.add_index(
                field_name=name,
                index_name=f"{name}_index",
                index_type="DISKANN",
                metric_type="COSINE",
            )
            schema.add_field(
                field_name=f"has_{name}",
                datatype=DataType.BOOL,
            )

        if include_sequence:
            schema.add_field(
                field_name="sequence",
                datatype=DataType.VARCHAR,
                max_length=65535,
            )
            schema.add_field(
                field_name="seq_length",
                datatype=DataType.INT32,
            )

        self.client.create_collection(
            collection_name=collection_name,
            schema=schema,
            index_params=index_params,
            num_shards=8,
        )

    def _initialize_collection_from_batch(
        self,
        collection_name: str,
        batch: EmbeddingBatch,
        include_sequence: bool,
    ) -> None:
        """Initialize collection schema from first EmbeddingBatch.

        This is a blocking operation that creates the collection if it doesn't exist.
        Uses the first batch to determine vector field names, dtypes, and dimensions.

        Args:
            collection_name: Name of the collection to create.
            batch: First EmbeddingBatch to use for schema inference.
            include_sequence: Whether to include sequence field in schema.
        """
        # Check if collection already exists
        if collection_name in self._initialized_collections:
            logger.warning(
                f"Collection '{collection_name}' already initialized. Skipping initialization"
            )
            return

        if not batch.embeddings:
            raise ValueError("Cannot initialize collection: batch has no embeddings")

        # Extract vector field information from batch
        vec_field_names: list[str] = []
        vec_field_dtypes: list[DataType] = []
        vec_dims: list[int] = []

        for pooling_name, embeddings in batch.embeddings.items():
            if not embeddings:
                continue

            # Get first embedding to determine dtype and dimension
            first_emb = embeddings[0]
            if not isinstance(first_emb, np.ndarray):
                raise ValueError(f"Expected numpy array for embedding, got {type(first_emb)}")

            # Determine dtype
            if first_emb.dtype == np.float16:
                dtype = DataType.FLOAT16_VECTOR
            elif first_emb.dtype == np.float32:
                dtype = DataType.FLOAT_VECTOR
            else:
                raise ValueError(f"Unsupported embedding dtype: {first_emb.dtype}")

            # Get dimension
            dim = int(first_emb.shape[0])

            vec_field_names.append(pooling_name)
            vec_field_dtypes.append(dtype)
            vec_dims.append(dim)

        if not vec_field_names:
            raise ValueError("Cannot initialize collection: no valid embeddings found in batch")

        logger.debug(
            f"Initializing collection '{collection_name}' with fields: {vec_field_names}, "
            f"dtypes: {vec_field_dtypes}, dimensions: {vec_dims}"
        )

        # Create collection
        self.create_collection(
            collection_name=collection_name,
            vec_field_names=vec_field_names,
            vec_field_dtypes=vec_field_dtypes,
            vec_dims=vec_dims,
            include_sequence=include_sequence,
        )

        self._initialized_collections.add(collection_name)
        logger.info(f"Collection '{collection_name}' initialized successfully")

    def _check_existing_ids(self, ids: list[str]) -> set[str]:
        """Check if the ids are already in the collection."""
        response = self.client.get(self.collection_name, ids=ids, output_fields=["protein_id"])
        return set(result["protein_id"] for result in response)

    async def insert_async(
        self,
        collection_name: str,
        batch: EmbeddingBatch,
        include_sequence: bool,
    ) -> int:
        """Insert EmbeddingBatch into Milvus collection asynchronously.

        Automatically initializes collection if it doesn't exist.
        Converts EmbeddingBatch to Milvus format and inserts using async client.

        Args:
            collection_name: Name of the collection to insert into.
            batch: EmbeddingBatch containing protein_ids, sequences, and embeddings.
            include_sequence: Whether to include sequence in inserted records.

        Returns:
            Number of records inserted.

        Raises:
            ValueError: If batch is empty or has mismatched lengths.
        """
        if not batch.protein_ids:
            logger.warning("Empty batch, skipping insert")
            return 0

        if len(batch.protein_ids) != len(batch.sequences):
            raise ValueError(
                f"Mismatched lengths: {len(batch.protein_ids)} protein_ids vs "
                f"{len(batch.sequences)} sequences"
            )

        # Initialize collection if needed (blocking, but only once)
        if collection_name not in self._initialized_collections:
            self._initialize_collection_from_batch(collection_name, batch, include_sequence)

        # Convert EmbeddingBatch to Milvus format
        milvus_data: list[dict[str, object]] = []

        for i, protein_id in enumerate(batch.protein_ids):
            record: dict[str, object] = {
                "protein_id": protein_id,
            }

            if include_sequence:
                sequence = batch.sequences[i]
                record["sequence"] = sequence
                record["seq_length"] = len(sequence)

            # Add embeddings for each pooling method
            for pooling_name, embeddings in batch.embeddings.items():
                if i >= len(embeddings):
                    logger.warning(f"Missing embedding for {pooling_name} at index {i}, skipping")
                    continue

                embedding = embeddings[i]
                if not isinstance(embedding, np.ndarray):
                    logger.warning(
                        f"Invalid embedding type for {pooling_name} at index {i}, skipping"
                    )
                    continue

                record[pooling_name] = embedding
                record[f"has_{pooling_name}"] = bool(embedding.any())

            milvus_data.append(record)

        if not milvus_data:
            logger.warning("No valid data to insert after conversion")
            return 0

        # Insert using sync client wrapped in thread pool
        # This avoids event loop issues with AsyncMilvusClient
        def _insert() -> None:
            self.client.insert(collection_name, milvus_data)

        await asyncio.to_thread(_insert)

        logger.debug(f"Inserted {len(milvus_data)} records into '{collection_name}'")
        return len(milvus_data)

    # Search
    # ------------------------------------------------------------

    def get_vectors(
        self, collection_name: str, protein_ids: list[str], vector_field_name: str = "mean_pooling"
    ) -> list[np.ndarray]:
        response = self.client.get(
            collection_name, ids=protein_ids, output_fields=[vector_field_name]
        )

        np_array = np.vstack([result[vector_field_name] for result in response]).astype(
            self._get_numpy_type(collection_name, vector_field_name)
        )
        return np_array

    def vector_search(
        self,
        collection_name: str,
        query_vectors: list[np.ndarray],
        vector_field_name: str = "mean_pooling",
        return_vector: bool = False,
        n_hits: int = 10,
        return_fields: list[str] = [],
        radius: float = 1.0,
        range_filter: float = 0.0,
        hitlist_size: int = 16384,
        offset: int = 0,
    ) -> list[list[dict[str, object]]]:
        if return_vector:
            return_fields.append(vector_field_name)

        return_fields = list(set(return_fields))

        # if vector is 1d make it 2d
        if len(query_vectors.shape) == 1:
            query_vectors = [query_vectors]

        records: list[dict]

        records = []
        records.extend(
            self.client.search(
                collection_name=collection_name,
                data=query_vectors,
                anns_field=vector_field_name,
                limit=n_hits,
                output_fields=return_fields,
                params={
                    "offset": offset,
                    "radius": radius,
                    "range_filter": range_filter,
                },
            )
        )

        # flatten the records
        clean_recs = []
        for rec in records[0]:
            value = rec["entity"][vector_field_name]
            vector = {
                vector_field_name: self._get_numpy_type(collection_name, vector_field_name)(value)
            }
            del rec["entity"]
            rec.update(vector)
            clean_recs.append(rec)

        df = pd.DataFrame.from_records(clean_recs)
        return df

    def id_search(self, collection_name: str, protein_ids: list[str]) -> list[dict[str, object]]:
        if isinstance(protein_ids, str):
            protein_ids = [protein_ids]

        # get

    def _get_search_load_params(
        self,
        n_hits: int,
        n_queries: int | None = None,
        *,
        max_total_hits_per_request: int = 16384,
        request_size_limit_mb: float = 64.0,
        safety: float = 0.6,
        id_bytes: int = 8,
        dist_bytes: int = 4,
        wire_overhead: float = 1.6,
    ) -> tuple[int, int]:
        """Return (hitlist_size, queries_per_search) without exceeding limits.

        n_hits: desired results per query (topK).
        n_queries: optional cap on queries per network request.
        """

        # Byte-budget-derived cap on total results we can return in one request.
        request_budget_bytes = int(request_size_limit_mb * (1024**2) * safety)
        bytes_per_result = int((id_bytes + dist_bytes) * wire_overhead)
        budget_cap_total_hits = request_budget_bytes // bytes_per_result

        # Effective total-results cap per request respects both byte budget and server max.
        total_hits_cap = min(budget_cap_total_hits, max_total_hits_per_request)

        # Per-query result count cannot exceed requested n_hits nor the total cap.
        hitlist_size = min(n_hits, total_hits_cap)

        # Number of queries we can pack while staying within the total results cap.
        queries_per_search = total_hits_cap // hitlist_size

        if n_queries is not None:
            queries_per_search = min(queries_per_search, n_queries)

        # Ensure product stays within the cap (guard against rounding).
        if queries_per_search * hitlist_size > total_hits_cap:
            queries_per_search = total_hits_cap // hitlist_size

        return hitlist_size, queries_per_search

    def get_all_vectors(
        self,
        collection_name: str,
        vector_field_name: str = "mean_pooling",
        show_progress: bool = True,
    ) -> list[np.ndarray]:
        it = self.client.query_iterator(
            collection_name=collection_name,
            batch_size=1000,
            output_fields=[vector_field_name],
        )
        protein_ids = []
        vectors = []

        # get total number of rows
        total_rows = self.client.get_collection_stats(collection_name)["row_count"]

        with create_progress() as progress:
            task = progress.add_task(
                f"Loading {vector_field_name} vectors from {collection_name}", total=total_rows
            )

            while True:
                batch = it.next()
                if not batch:
                    break

                # extract protein id and make in list of protein ids
                for rec in batch:
                    protein_ids.append(rec["protein_id"])
                    vectors.append(rec[vector_field_name])

                progress.update(task, advance=len(batch))

            task_convert = progress.add_task(
                "Converting vectors to numpy array", total=None, start=True
            )
            vectors = np.vstack(vectors).astype(
                self._get_numpy_type(collection_name, vector_field_name)
            )
            progress.update(task_convert, completed=True)
            sleep(1)

            return protein_ids, vectors

    def _get_numpy_type(
        self,
        collection_name: str,
        vector_field_name: str,
    ) -> np.dtype:
        schema = self.client.describe_collection(collection_name)
        for field in schema["fields"]:
            if field["name"] == vector_field_name:
                return MILVUS_TO_NP_DTYPE_MAP[field["type"]]
        raise ValueError(
            f"Vector field '{vector_field_name}' not found in collection '{collection_name}'"
        )
