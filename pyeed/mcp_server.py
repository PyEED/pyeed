from __future__ import annotations

from fastmcp import FastMCP
from toon_format import encode

from pyeed.db.milvus import get_async_milvus_client
from pyeed.db.neo4j import get_async_driver
from pyeed.ingest.model.protein import Protein
from pyeed.queries import get_similar_proteins_by_ids

mcp = FastMCP("Pyeed MCP Server")

driver = get_async_driver()
milvus_client = get_async_milvus_client()


@mcp.tool
async def get_proteins(ids: list[str]) -> str:
    """Get proteins from the database by IDs."""

    proteins = await Protein.get_filtered(driver, ids=ids)
    return encode([p.model_dump(exclude_unset=True) for p in proteins])


@mcp.tool
async def protein_similarity_search(ids: list[str], limit: int = 10) -> str:
    """Search for similar proteins by IDs.

    Args:
        ids: List of protein IDs to search for.
        limit: Maximum number of similar proteins to return (default: 10).

    """
    proteins = await get_similar_proteins_by_ids(
        ids=ids,
        neo4j_driver=driver,
        milvus_client=milvus_client,
        collection_name="pyeed",
        limit=limit,
        vector_field_name="mean_pooling",
    )
    return encode([p.model_dump(exclude_unset=True) for p in proteins])


if __name__ == "__main__":
    mcp.run(transport="sse")
