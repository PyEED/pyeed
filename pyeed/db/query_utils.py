from __future__ import annotations

from collections.abc import Awaitable, Callable
from typing import Any

from neo4j import AsyncDriver, AsyncManagedTransaction

# --------- Result Processors --------- #


async def process_single_record(
    tx: AsyncManagedTransaction,
    query: str,
    params: dict[str, Any],
) -> dict[str, Any] | None:
    """Process a query that returns a single record."""
    result = await tx.run(query, params)
    record = await result.single()
    return dict(record) if record else None


async def process_multiple_records(
    tx: AsyncManagedTransaction,
    query: str,
    params: dict[str, Any],
) -> list[dict[str, Any]]:
    """Process a query that returns multiple records."""
    result = await tx.run(query, params)
    return [dict(record) async for record in result]


async def process_count(
    tx: AsyncManagedTransaction,
    query: str,
    params: dict[str, Any],
) -> int:
    """Process a count query."""
    result = await tx.run(query, params)
    record = await result.single()
    return record[0] if record else 0


# --------- Transaction Executors --------- #

type ResultProcessor[T] = Callable[[AsyncManagedTransaction, str, dict[str, Any]], Awaitable[T]]


async def execute_read[T](
    driver: AsyncDriver,
    query: str,
    params: dict[str, Any],
    processor: ResultProcessor[T],
) -> T:
    """Execute a read transaction with proper scope handling.

    - Session is properly closed via context manager
    - Results are consumed within transaction scope
    - Generic return type via processor function
    """
    async with driver.session() as session:
        return await session.execute_read(processor, query, params)


async def execute_write(
    driver: AsyncDriver,
    query: str,
    params: dict[str, Any] | None = None,
    rows: list[dict[str, Any]] | None = None,
    session_kwargs: dict[str, Any] | None = None,
) -> None:
    """Execute a write transaction.

    Supports both simple parameterized queries and UNWIND batch operations.

    Args:
        driver: Neo4j async driver.
        query: Cypher query string.
        params: Simple parameter dictionary (for regular queries).
        rows: List of dictionaries for UNWIND batch operations.
        session_kwargs: Optional session configuration.

    Examples:
        # Simple write
        await execute_write(driver, "CREATE (n:Node {id: $id})", params={"id": "123"})

        # Batch write with UNWIND
        await execute_write(
            driver,
            "UNWIND $rows AS row MERGE (n:Node {id: row.id})",
            rows=[{"id": "1"}, {"id": "2"}]
        )
    """
    if params is not None and rows is not None:
        raise ValueError("Provide either 'params' or 'rows', not both")

    async def _tx(
        tx: AsyncManagedTransaction, q: str, p: dict[str, Any] | list[dict[str, Any]]
    ) -> None:
        await tx.run(q, p)

    kwargs = session_kwargs or {}
    async with driver.session(**kwargs) as session:
        payload = rows if rows is not None else (params or {})
        await session.execute_write(_tx, query, payload)
