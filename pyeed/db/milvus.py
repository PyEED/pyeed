from dataclasses import dataclass

from loguru import logger
from pymilvus import (
    AsyncMilvusClient,
    Collection,
    DataType,
    MilvusClient,
)


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
            max_length=64,
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
    vector_db.client.drop_collection("protein_emb")

    # create collection
    print("Creating collection 'protein_emb' with vector fields...")
    vector_db.create_collection(
        collection_name="protein_emb",
        vec_field_names=["vec_mean_pooling"],
        vec_field_dtypes=[DataType.FLOAT16_VECTOR],
        vec_dims=[1024],
        include_sequence=True,
    )

    # check collection
    print("Listing collections after creation of 'protein_emb':")
    print(f"Collections: {vector_db.client.list_collections()}")
    print("Describing 'protein_emb' collection schema:")
    print(f"Collection schema: {vector_db.client.describe_collection('protein_emb')}")

    import numpy as np
    from rich import print

    # insert data
    print("Populating data with 100 protein entries and adding one zeroed entry...")
    data = []
    for i in range(100):
        data.append(
            {
                "protein_id": f"P0000{i}",
                "vec_mean_pooling": np.random.RandomState(42).rand(1024).astype(np.float16),
                "sequence": "ACDEFGHIKLMNPQRSTVWY",
                "seq_length": len("ACDEFGHIKLMNPQRSTVWY"),
                "has_vec_mean_pooling": True,
            }
        )

    # add entry where embedding and sequece is none
    print("Adding one entry with zeroed embedding and short sequence.")
    data.append(
        {
            "protein_id": "P69",
            "vec_mean_pooling": np.zeros(1024, dtype=np.float16),
            "sequence": "ASD",
            "seq_length": 3,
            "has_vec_mean_pooling": False,
        }
    )
    print(f"Upserting {len(data)} data entries into 'protein_emb' collection...")
    vector_db.client.insert("protein_emb", data)
    print("Upsert finished.")
    vector_db.client.flush("protein_emb")  # force persistence

    # first vec
    first_vec = data[0]["vec_mean_pooling"]
    print("The dtype of the first vector inserted:")
    print(first_vec.dtype)

    # get embeddig and protein id of P00000
    print('Querying for protein with id "P00000"...')
    res = vector_db.client.query(
        collection_name="protein_emb",
        filter='protein_id == "P69"',
        output_fields=[
            "protein_id",
            "seq_length",
            "sequence",
            "vec_mean_pooling",
            "has_vec_mean_pooling",
        ],  # include vector if you want
    )
    print(f"Number of results for protein_id == 'P00000': {len(res)}")
    row = res[0]  # one entity from client.query(...) or client.search(...)
    print("Fields in the first result row:", row.keys())
    print("First result row:")
    print(row)
    buf = row["vec_mean_pooling"][0]  # this is a bytes-like object
    print("Converting 'vec_mean_pooling' from buffer to numpy float16 array...")
    vec_f16 = np.frombuffer(buf, dtype=np.float16)
    print("Retrieved vector (float16):")
    print(vec_f16)
    print("Original (first) vector inserted:")
    print(first_vec)
    print("Are the retrieved vector and first vector equal?")
    print(vec_f16 == first_vec)
