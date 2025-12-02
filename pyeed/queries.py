from collections.abc import Awaitable, Callable
from typing import Any, Literal

import numpy as np
import numpy.typing as npt
from neo4j import AsyncDriver, AsyncManagedTransaction
from pydantic import BaseModel
from pymilvus import AsyncMilvusClient, DataType

from pyeed.db.milvus import get_async_milvus_client
from pyeed.db.neo4j import get_async_driver
from pyeed.ingest.model import Protein

MILVUS_TO_NP_DTYPE_MAP = {
    DataType.FLOAT16_VECTOR: np.float16,
    DataType.FLOAT_VECTOR: np.float32,
}


class VectorResponse(BaseModel):
    id: str
    vector: list[float]
    dtype: Literal["float16", "float32"] = "float32"
    dim: int


# Generic result processors
async def process_single_record(
    tx: AsyncManagedTransaction, query: str, params: dict[str, Any]
) -> dict[str, Any] | None:
    """Process a query that returns a single record."""
    result = await tx.run(query, params)
    record = await result.single()
    return dict(record) if record else None


async def process_multiple_records(
    tx: AsyncManagedTransaction, query: str, params: dict[str, Any]
) -> list[dict[str, Any]]:
    """Process a query that returns multiple records."""
    result = await tx.run(query, params)
    return [dict(record) async for record in result]


async def process_count(tx: AsyncManagedTransaction, query: str, params: dict[str, Any]) -> int:
    """Process a count query."""
    result = await tx.run(query, params)
    record = await result.single()
    return record[0] if record else 0


async def execute_read_transaction(
    neo4j_driver: AsyncDriver,
    query: str,
    params: dict[str, Any],
    result_processor: Callable[[AsyncManagedTransaction, str, dict[str, Any]], Awaitable[Any]],
) -> Any:
    """Generic method to execute read transactions with custom result processing."""
    async with neo4j_driver.session() as session:
        return await session.execute_read(result_processor, query, params)


async def get_protein_by_id(
    id: str,
    neo4j_driver: AsyncDriver,
) -> Protein | None:
    query = "MATCH (p:Protein {id: $id}) RETURN properties(p) AS protein"

    protein_data = await execute_read_transaction(
        neo4j_driver=neo4j_driver,
        query=query,
        params={"id": id},
        result_processor=process_single_record,
    )

    if protein_data and "protein" in protein_data:
        return Protein(**protein_data["protein"])
    return None


async def get_all_proteins(neo4j_driver: AsyncDriver) -> list[Protein]:
    query = "MATCH (p:Protein) RETURN properties(p) AS protein"

    records = await execute_read_transaction(
        neo4j_driver=neo4j_driver,
        query=query,
        params={},
        result_processor=process_multiple_records,
    )

    return [Protein(**record["protein"]) for record in records]


async def count_nodes_per_label(neo4j_driver: AsyncDriver) -> dict[str, int]:
    """Count the number of nodes for each label."""
    query = "MATCH (n) RETURN labels(n) AS labels, count(n) AS count"

    records = await execute_read_transaction(
        neo4j_driver=neo4j_driver,
        query=query,
        params={},
        result_processor=process_multiple_records,
    )

    return {entry["labels"][0]: entry["count"] for entry in records}


# ------------------- Similarity Search -------------------
async def get_vectors_by_ids(
    ids: list[str],
    milvus_client: AsyncMilvusClient,
    collection_name: str,
    vector_field_name: str = "mean_pooling",
) -> dict[str, npt.NDArray[np.float32] | npt.NDArray[np.float16]]:
    """Get vectors by IDs.

    Returns:
        A dictionary of protein IDs and their vectors.
    """
    records = await milvus_client.get(
        collection_name=collection_name,
        ids=ids,
        output_fields=[vector_field_name],
    )
    schema = await milvus_client.describe_collection(collection_name)
    for field in schema["fields"]:
        if field["name"] == vector_field_name:
            dtype = MILVUS_TO_NP_DTYPE_MAP[field["type"]]
            break
    else:
        raise ValueError(
            f"Vector field '{vector_field_name}' not found in collection '{collection_name}'"
        )

    return {
        record["protein_id"]: np.array(record[vector_field_name], dtype=dtype) for record in records
    }


class MilvusSearchResult(BaseModel):
    query_id: str
    target_id: str
    distance: float
    vector: list[float]
    dtype: Literal["float16", "float32"]


async def get_similar_proteins_by_ids(
    ids: list[str],
    neo4j_driver: AsyncDriver,
    milvus_client: AsyncMilvusClient,
    collection_name: str,
    limit: int,
    vector_field_name: str,
) -> list[Protein]:
    """Get similar proteins by ID."""

    if isinstance(ids, str):
        ids = [ids]

    vectors = await get_vectors_by_ids(
        ids=ids,
        milvus_client=milvus_client,
        collection_name=collection_name,
        vector_field_name="mean_pooling",
    )

    if not vectors:
        return []

    np_dtype = vectors[next(iter(vectors.keys()))].dtype.name

    records = await milvus_client.search(
        collection_name=collection_name,
        data=list(vectors.values()),
        limit=limit,
        anns_field=vector_field_name,
        output_fields=["protein_id", vector_field_name],
    )

    search_results = []
    for record_idx, record in enumerate(records):
        for result in record:
            search_results.append(
                MilvusSearchResult(
                    query_id=ids[record_idx],
                    target_id=result["protein_id"],
                    distance=result["distance"],
                    vector=result[vector_field_name],
                    dtype=np_dtype,
                )
            )

    return search_results


async def main() -> None:
    from rich import print

    neo4j_driver = get_async_driver()
    milvus_client = get_async_milvus_client()
    res = await get_similar_proteins_by_ids(
        ids=["P12345", "P0CW62"],
        neo4j_driver=neo4j_driver,
        milvus_client=milvus_client,
        collection_name="pyeed",
        limit=10,
        vector_field_name="mean_pooling",
    )
    print(res)

    await neo4j_driver.close()
    await milvus_client.close()

    print("Done")


if __name__ == "__main__":
    import asyncio

    asyncio.run(main())
