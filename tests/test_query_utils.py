"""Unit tests for Neo4j query utilities."""

from unittest.mock import AsyncMock, MagicMock

import pytest

from pyeed.db.query_utils import (
    execute_read,
    execute_write,
    process_count,
    process_multiple_records,
    process_single_record,
)

# ============================================================================
# Test Fixtures
# ============================================================================


@pytest.fixture
def mock_neo4j_transaction():
    """Create a mock Neo4j async transaction."""
    tx = AsyncMock()
    result = AsyncMock()
    tx.run = AsyncMock(return_value=result)
    return tx, result


@pytest.fixture
def mock_neo4j_driver():
    """Create a mock Neo4j async driver."""
    driver = MagicMock()
    session = AsyncMock()

    # Create async context manager for session
    async_context_manager = MagicMock()
    async_context_manager.__aenter__ = AsyncMock(return_value=session)
    async_context_manager.__aexit__ = AsyncMock(return_value=None)

    driver.session = MagicMock(return_value=async_context_manager)

    return driver, session


# ============================================================================
# Test Result Processors
# ============================================================================


class TestProcessSingleRecord:
    """Tests for process_single_record processor."""

    @pytest.mark.asyncio
    async def test_process_single_record_with_data(self, mock_neo4j_transaction):
        """Test processing a single record that exists."""
        tx, result = mock_neo4j_transaction

        # Mock record with data - Neo4j records are dict-like
        # When dict(record) is called, it iterates over key-value pairs
        mock_record = {"key1": "value1", "key2": "value2"}
        result.single = AsyncMock(return_value=mock_record)

        query = "MATCH (n) RETURN n"
        params = {}

        result_data = await process_single_record(tx, query, params)

        assert result_data == {"key1": "value1", "key2": "value2"}
        tx.run.assert_called_once_with(query, params)

    @pytest.mark.asyncio
    async def test_process_single_record_no_data(self, mock_neo4j_transaction):
        """Test processing when no record exists."""
        tx, result = mock_neo4j_transaction

        result.single = AsyncMock(return_value=None)

        query = "MATCH (n) RETURN n"
        params = {}

        result_data = await process_single_record(tx, query, params)

        assert result_data is None


class TestProcessMultipleRecords:
    """Tests for process_multiple_records processor."""

    @pytest.mark.asyncio
    async def test_process_multiple_records(self, mock_neo4j_transaction):
        """Test processing multiple records."""
        tx, result = mock_neo4j_transaction

        # Mock multiple records - Neo4j records are dict-like
        mock_record1 = {"id": "1", "name": "One"}
        mock_record2 = {"id": "2", "name": "Two"}

        async def async_iter():
            yield mock_record1
            yield mock_record2

        result.__aiter__ = lambda self: async_iter()

        query = "MATCH (n) RETURN n"
        params = {}

        result_data = await process_multiple_records(tx, query, params)

        assert len(result_data) == 2
        assert result_data[0] == {"id": "1", "name": "One"}
        assert result_data[1] == {"id": "2", "name": "Two"}
        tx.run.assert_called_once_with(query, params)

    @pytest.mark.asyncio
    async def test_process_multiple_records_empty(self, mock_neo4j_transaction):
        """Test processing when no records exist."""
        tx, result = mock_neo4j_transaction

        async def async_iter():
            return
            yield  # Make it an async generator

        result.__aiter__ = lambda self: async_iter()

        query = "MATCH (n) RETURN n"
        params = {}

        result_data = await process_multiple_records(tx, query, params)

        assert result_data == []


class TestProcessCount:
    """Tests for process_count processor."""

    @pytest.mark.asyncio
    async def test_process_count_with_value(self, mock_neo4j_transaction):
        """Test processing a count query with result."""
        tx, result = mock_neo4j_transaction

        mock_record = MagicMock()
        mock_record.__getitem__ = lambda self, key: 42 if key == 0 else None
        result.single = AsyncMock(return_value=mock_record)

        query = "MATCH (n) RETURN count(n)"
        params = {}

        count = await process_count(tx, query, params)

        assert count == 42
        tx.run.assert_called_once_with(query, params)

    @pytest.mark.asyncio
    async def test_process_count_no_record(self, mock_neo4j_transaction):
        """Test processing count when no record exists."""
        tx, result = mock_neo4j_transaction

        result.single = AsyncMock(return_value=None)

        query = "MATCH (n) RETURN count(n)"
        params = {}

        count = await process_count(tx, query, params)

        assert count == 0


# ============================================================================
# Test Transaction Executors
# ============================================================================


class TestExecuteRead:
    """Tests for execute_read function."""

    @pytest.mark.asyncio
    async def test_execute_read_calls_session_execute_read(self, mock_neo4j_driver):
        """Test that execute_read properly calls session.execute_read."""
        driver, session = mock_neo4j_driver

        async def mock_processor(tx, query, params):
            return {"result": "data"}

        session.execute_read = AsyncMock(return_value={"result": "data"})

        query = "MATCH (n) RETURN n"
        params = {"param1": "value1"}

        result = await execute_read(driver, query, params, mock_processor)

        assert result == {"result": "data"}
        session.execute_read.assert_called_once()
        call_args = session.execute_read.call_args
        assert call_args[0][1] == query
        assert call_args[0][2] == params

    @pytest.mark.asyncio
    async def test_execute_read_creates_session(self, mock_neo4j_driver):
        """Test that execute_read creates and closes session."""
        driver, session = mock_neo4j_driver

        async def mock_processor(tx, query, params):
            return {"result": "data"}

        session.execute_read = AsyncMock(return_value={"result": "data"})

        await execute_read(driver, "MATCH (n) RETURN n", {}, mock_processor)

        # Verify session context manager was used
        driver.session.assert_called_once()


class TestExecuteWrite:
    """Tests for execute_write function."""

    @pytest.mark.asyncio
    async def test_execute_write_with_params(self, mock_neo4j_driver):
        """Test execute_write with params (simple query)."""
        driver, session = mock_neo4j_driver

        session.execute_write = AsyncMock()

        query = "CREATE (n:Node {id: $id})"
        params = {"id": "123"}

        await execute_write(driver, query, params=params)

        session.execute_write.assert_called_once()
        call_args = session.execute_write.call_args
        # Verify the transaction function was called
        assert call_args[0][1] == query
        assert call_args[0][2] == params

    @pytest.mark.asyncio
    async def test_execute_write_with_rows(self, mock_neo4j_driver):
        """Test execute_write with rows (UNWIND batch operation)."""
        driver, session = mock_neo4j_driver

        session.execute_write = AsyncMock()

        query = "UNWIND $rows AS row MERGE (n:Node {id: row.id})"
        rows = [{"id": "1"}, {"id": "2"}]

        await execute_write(driver, query, rows=rows)

        session.execute_write.assert_called_once()
        call_args = session.execute_write.call_args
        assert call_args[0][1] == query
        assert call_args[0][2] == rows

    @pytest.mark.asyncio
    async def test_execute_write_with_session_kwargs(self, mock_neo4j_driver):
        """Test execute_write passes session_kwargs to session."""
        driver, session = mock_neo4j_driver

        session.execute_write = AsyncMock()

        query = "CREATE (n:Node)"
        session_kwargs = {"database": "test_db"}

        await execute_write(driver, query, params={}, session_kwargs=session_kwargs)

        # Verify session was created with kwargs
        driver.session.assert_called_once_with(**session_kwargs)

    @pytest.mark.asyncio
    async def test_execute_write_raises_error_both_params_and_rows(self, mock_neo4j_driver):
        """Test execute_write raises error when both params and rows provided."""
        driver, session = mock_neo4j_driver

        with pytest.raises(ValueError, match="Provide either 'params' or 'rows', not both"):
            await execute_write(driver, "CREATE (n:Node)", params={"id": "1"}, rows=[{"id": "1"}])

    @pytest.mark.asyncio
    async def test_execute_write_with_empty_params(self, mock_neo4j_driver):
        """Test execute_write with empty params dict."""
        driver, session = mock_neo4j_driver

        session.execute_write = AsyncMock()

        query = "MATCH (n) RETURN n"

        await execute_write(driver, query, params={})

        session.execute_write.assert_called_once()
        call_args = session.execute_write.call_args
        assert call_args[0][2] == {}

    @pytest.mark.asyncio
    async def test_execute_write_with_none_params(self, mock_neo4j_driver):
        """Test execute_write with None params (uses empty dict)."""
        driver, session = mock_neo4j_driver

        session.execute_write = AsyncMock()

        query = "MATCH (n) RETURN n"

        await execute_write(driver, query, params=None)

        session.execute_write.assert_called_once()
        call_args = session.execute_write.call_args
        assert call_args[0][2] == {}

    @pytest.mark.asyncio
    async def test_execute_write_transaction_runs_query(self, mock_neo4j_driver):
        """Test that the internal transaction function runs the query."""
        driver, session = mock_neo4j_driver

        mock_tx = AsyncMock()
        captured_query = None
        captured_payload = None

        async def mock_execute_write(tx_func, query, payload):
            nonlocal captured_query, captured_payload
            captured_query = query
            captured_payload = payload
            await tx_func(mock_tx, query, payload)

        session.execute_write = mock_execute_write

        query = "CREATE (n:Node {id: $id})"
        params = {"id": "123"}

        await execute_write(driver, query, params=params)

        assert captured_query == query
        assert captured_payload == params
        mock_tx.run.assert_called_once_with(query, params)
