from __future__ import annotations

import asyncio
from dataclasses import dataclass

import numpy as np
from loguru import logger
from pymilvus import (
    AsyncMilvusClient,
    Collection,
    DataType,
    MilvusClient,
)

from ..embedding.types import EmbeddingBatch


@dataclass
class UploadStats:
    """Statistics for upload operations."""

    total_inserted: int = 0
    total_batches: int = 0
    failed_inserts: int = 0


class VectorDB:
    def __init__(
        self,
        uri: str | None = None,
        token: str | None = None,
        batch_size: int = 1000,
        max_batch_mb: float = 30.0,
    ):
        self.batch_size = batch_size
        self.max_batch_mb = max_batch_mb

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
            enable_dynamic_field=True,
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
        include_sequence: bool = True,
    ) -> None:
        """Initialize collection schema from first EmbeddingBatch.

        This is a blocking operation that creates the collection if it doesn't exist.
        Uses the first batch to determine vector field names, dtypes, and dimensions.

        Args:
            collection_name: Name of the collection to create.
            batch: First EmbeddingBatch to use for schema inference.
            include_sequence: Whether to include sequence field in schema.
        """
        if collection_name in self._initialized_collections:
            return

        # Check if collection already exists
        if collection_name in self.client.list_collections():
            self._initialized_collections.add(collection_name)
            logger.info(f"Collection '{collection_name}' already exists, skipping initialization")
            return

        if not batch.embeddings:
            raise ValueError("Cannot initialize collection: batch has no embeddings")

        # Extract vector field information from batch
        vec_field_names: list[str] = []
        vec_field_dtypes: list[DataType] = []
        vec_dims: list[int] = []

        for pool_name, embeddings in batch.embeddings.items():
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

            # Field name format: vec_{pool_name}
            field_name = f"vec_{pool_name}"
            vec_field_names.append(field_name)
            vec_field_dtypes.append(dtype)
            vec_dims.append(dim)

        if not vec_field_names:
            raise ValueError("Cannot initialize collection: no valid embeddings found in batch")

        logger.info(
            f"Initializing collection '{collection_name}' with {len(vec_field_names)} vector fields"
        )

        # Create collection (blocking operation)
        self.create_collection(
            collection_name=collection_name,
            vec_field_names=vec_field_names,
            vec_field_dtypes=vec_field_dtypes,
            vec_dims=vec_dims,
            include_sequence=include_sequence,
        )

        self._initialized_collections.add(collection_name)
        logger.info(f"Collection '{collection_name}' initialized successfully")

    async def insert_async(
        self,
        collection_name: str,
        batch: EmbeddingBatch,
        include_sequence: bool = True,
    ) -> int:
        """Insert EmbeddingBatch into Milvus collection asynchronously.

        Automatically initializes collection if it doesn't exist (blocking operation).
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
            for pool_name, embeddings in batch.embeddings.items():
                if i >= len(embeddings):
                    logger.warning(f"Missing embedding for {pool_name} at index {i}, skipping")
                    continue

                embedding = embeddings[i]
                if not isinstance(embedding, np.ndarray):
                    logger.warning(f"Invalid embedding type for {pool_name} at index {i}, skipping")
                    continue

                field_name = f"vec_{pool_name}"
                record[field_name] = embedding
                record[f"has_{field_name}"] = True

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


# Example usage
if __name__ == "__main__":
    from rich import print

    print("Initializing VectorDB...")
    vector_db = VectorDB(
        uri="http://localhost:19530",
        token="root:Milvus",
    )

    # check connection
    print("Checking available collections and databases after connection.")
    print(f"Collections: {vector_db.collections}")
    print(f"Databases: {vector_db.databases}")

    # drop collection
    print("Dropping collection 'protein_emb' if it exists...")
    vector_db.client.drop_collection("test")
