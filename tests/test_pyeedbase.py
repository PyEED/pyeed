"""Unit tests for BaseNode and LabelProperty validation."""

from typing import Annotated
from unittest.mock import AsyncMock, MagicMock, patch

import pytest
from pydantic import Field

from pyeed.ingest.model.pyeedbase import BaseNode, LabelProperty

# ============================================================================
# Test Fixtures
# ============================================================================


@pytest.fixture
def valid_node_class():
    """Create a valid node class with single indexed field."""

    class TestNode(BaseNode):
        id: Annotated[str, LabelProperty(index=True)] = Field(..., description="ID")
        name: str = Field(..., description="Name")
        age: int | None = Field(None, description="Age")

    return TestNode


@pytest.fixture
def sample_node(valid_node_class):
    """Create a sample node instance."""
    return valid_node_class(id="test-123", name="Test Node", age=42)


# ============================================================================
# Test LabelProperty Index Validation
# ============================================================================


class TestLabelPropertyValidation:
    """Tests for LabelProperty index validation at class definition time."""

    def test_single_indexed_field_allowed(self, valid_node_class):
        """Test that a class with single indexed field is allowed."""
        # Should not raise
        node = valid_node_class(id="test", name="Node")
        assert node.id == "test"
        assert node.name == "Node"

    def test_multiple_indexed_fields_raises_error(self):
        """Test that multiple indexed fields raise ValueError at class definition."""
        with pytest.raises(ValueError, match="Only one field can be marked"):

            class MultiIndexNode(BaseNode):
                id: Annotated[str, LabelProperty(index=True)] = Field(..., description="ID")
                code: Annotated[str, LabelProperty(index=True)] = Field(..., description="Code")

            MultiIndexNode(id="test", code="123")

    def test_no_indexed_field_allowed_at_creation(self):
        """Test that class without indexed field can be created."""

        class NoIndexNode(BaseNode):
            name: str = Field(..., description="Name")

        # Should create successfully
        node = NoIndexNode(name="Test")
        assert node.name == "Test"

    def test_error_message_includes_field_names(self):
        """Test that error message includes field names."""
        with pytest.raises(ValueError, match=r"id.*code"):

            class BadNode(BaseNode):
                id: Annotated[str, LabelProperty(index=True)] = Field(..., description="ID")
                code: Annotated[str, LabelProperty(index=True)] = Field(..., description="Code")

            BadNode(id="test", code="123")


# ============================================================================
# Test get_index_field Method
# ============================================================================


class TestGetIndexField:
    """Tests for get_index_field() method."""

    def test_returns_indexed_field_name(self, valid_node_class):
        """Test that get_index_field returns the indexed field name."""
        assert valid_node_class.get_index_field() == "id"

    def test_instance_method_works(self, sample_node):
        """Test that get_index_field works on instance."""
        assert sample_node.get_index_field() == "id"

    def test_no_indexed_field_raises_error(self):
        """Test that class without indexed field raises ValueError."""

        class NoIndexNode(BaseNode):
            name: str = Field(..., description="Name")

        with pytest.raises(ValueError, match="No field with LabelProperty"):
            NoIndexNode.get_index_field()

    def test_error_message_helpful(self):
        """Test that error message provides guidance."""

        class NoIndexNode(BaseNode):
            name: str = Field(..., description="Name")

        with pytest.raises(
            ValueError, match="To mark a field as indexed, use: Annotated.*LabelProperty"
        ):
            NoIndexNode.get_index_field()


# ============================================================================
# Test Neo4j Integration (Mocked)
# ============================================================================


@pytest.fixture
def mock_neo4j_driver():
    """Create a mock Neo4j async driver."""
    driver = MagicMock()
    session = AsyncMock()

    # Create async context manager for session
    async_context_manager = MagicMock()
    async_context_manager.__aenter__ = AsyncMock(return_value=session)
    async_context_manager.__aexit__ = AsyncMock(return_value=None)

    # driver.session() returns the async context manager
    driver.session = MagicMock(return_value=async_context_manager)

    session.execute_write = AsyncMock()

    return driver


@pytest.fixture
def mock_neo4j_transaction():
    """Create a mock Neo4j async transaction."""
    tx = AsyncMock()
    tx.run = AsyncMock()
    return tx


class TestBulkUpsert:
    """Tests for _bulk_upsert() method with mocked Neo4j."""

    @pytest.mark.asyncio
    async def test__bulk_upsert_basic(self, valid_node_class, mock_neo4j_driver):
        """Test basic bulk upsert operation."""
        nodes = [
            valid_node_class(id="id1", name="Node 1"),
            valid_node_class(id="id2", name="Node 2"),
        ]

        await valid_node_class._bulk_upsert(mock_neo4j_driver, nodes)

        # Verify driver.session was called
        mock_neo4j_driver.session.assert_called_once()

    @pytest.mark.asyncio
    async def test_bulk_upsert_with_session_kwargs(self, valid_node_class, mock_neo4j_driver):
        """Test bulk upsert with session kwargs."""
        nodes = [valid_node_class(id="id1", name="Node 1")]
        session_kwargs = {"database": "test_db"}

        await valid_node_class._bulk_upsert(mock_neo4j_driver, nodes, session_kwargs=session_kwargs)

        # Verify session was called with kwargs
        mock_neo4j_driver.session.assert_called_once_with(**session_kwargs)

    @pytest.mark.asyncio
    async def test_bulk_upsert_empty_list(self, valid_node_class, mock_neo4j_driver):
        """Test bulk upsert with empty list does nothing."""
        await valid_node_class._bulk_upsert(mock_neo4j_driver, [])

        # Should not call session
        mock_neo4j_driver.session.assert_not_called()

    @pytest.mark.asyncio
    async def test_bulk_upsert_generates_correct_cypher(self, valid_node_class, mock_neo4j_driver):
        """Test that bulk upsert generates correct Cypher query."""
        nodes = [valid_node_class(id="id1", name="Node 1")]

        # Capture arguments passed to execute_write
        captured_query = None
        captured_rows = None

        async def mock_execute_write(driver, query, params=None, rows=None, session_kwargs=None):
            nonlocal captured_query, captured_rows
            captured_query = query
            captured_rows = rows

        with patch("pyeed.ingest.model.pyeedbase.execute_write", side_effect=mock_execute_write):
            await valid_node_class._bulk_upsert(mock_neo4j_driver, nodes)

        # Verify Cypher contains expected elements
        assert captured_query is not None
        assert "UNWIND $rows AS row" in captured_query
        assert "MERGE" in captured_query
        assert "TestNode" in captured_query
        assert "id: row.id" in captured_query
        assert "SET n += row" in captured_query

        # Verify rows data
        assert captured_rows is not None
        assert len(captured_rows) == 1
        assert captured_rows[0]["id"] == "id1"
        assert captured_rows[0]["name"] == "Node 1"


class TestUpsert:
    """Tests for single node upsert() method."""

    @pytest.mark.asyncio
    async def test_upsert_calls_bulk_upsert(self, sample_node, mock_neo4j_driver):
        """Test that upsert calls bulk_upsert with single node."""
        with patch.object(type(sample_node), "_bulk_upsert", new_callable=AsyncMock) as mock_bulk:
            await sample_node._upsert(mock_neo4j_driver)

            # Verify bulk_upsert was called with single node in list
            mock_bulk.assert_called_once()
            call_args = mock_bulk.call_args
            assert call_args[0][0] == mock_neo4j_driver
            assert len(call_args[0][1]) == 1
            assert call_args[0][1][0] == sample_node

    @pytest.mark.asyncio
    async def test_upsert_with_session_kwargs(self, sample_node, mock_neo4j_driver):
        """Test upsert passes session kwargs."""
        session_kwargs = {"database": "test_db"}

        with patch.object(type(sample_node), "_bulk_upsert", new_callable=AsyncMock) as mock_bulk:
            await sample_node._upsert(mock_neo4j_driver, session_kwargs=session_kwargs)

            # Verify session_kwargs was passed
            call_args = mock_bulk.call_args
            assert call_args[1]["session_kwargs"] == session_kwargs


# ============================================================================
# Integration Tests
# ============================================================================


class TestIntegration:
    """Integration tests combining multiple features."""

    def test_complete_workflow(self, valid_node_class):
        """Test complete workflow: create, validate, convert to dict."""
        # Create node
        node = valid_node_class(
            id="test-123",
            name="Integration Test",
            age=99,
        )

        # Verify node creation
        assert node.id == "test-123"
        assert node.name == "Integration Test"
        assert node.age == 99

        # Verify get_index_field
        assert node.get_index_field() == "id"

        # Verify model_dump includes everything
        d = node.model_dump()
        assert d["id"] == "test-123"
        assert d["name"] == "Integration Test"
        assert d["age"] == 99

    def test_inheritance_works(self):
        """Test that inheritance from BaseNode works correctly."""

        class ParentNode(BaseNode):
            id: Annotated[str, LabelProperty(index=True)] = Field(..., description="ID")
            parent_field: str = Field(..., description="Parent field")

        class ChildNode(ParentNode):
            child_field: str = Field(..., description="Child field")

        # Create child instance
        child = ChildNode(id="child-1", parent_field="parent", child_field="child")

        assert child.id == "child-1"
        assert child.parent_field == "parent"
        assert child.child_field == "child"

        # Verify get_index_field works on child
        assert child.get_index_field() == "id"

        # Verify model_dump includes all fields
        d = child.model_dump()
        assert d["id"] == "child-1"
        assert d["parent_field"] == "parent"
        assert d["child_field"] == "child"


# ============================================================================
# Test Query Methods
# ============================================================================


class TestQueryMethods:
    """Tests for query methods (get, get_by, get_all, count)."""

    @pytest.mark.asyncio
    async def test_get_single_id(self, valid_node_class, mock_neo4j_driver):
        """Test get() with single id parameter."""
        mock_data = {"node": {"id": "test-123", "name": "Test Node", "age": 42}}

        async def mock_execute_read(driver, query, params, processor):
            return mock_data

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            result = await valid_node_class.get(mock_neo4j_driver, id="test-123")

        assert result is not None
        assert result.id == "test-123"
        assert result.name == "Test Node"
        assert result.age == 42

    @pytest.mark.asyncio
    async def test_get_single_id_not_found(self, valid_node_class, mock_neo4j_driver):
        """Test get() returns None when node not found."""

        async def mock_execute_read(driver, query, params, processor):
            return None

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            result = await valid_node_class.get(mock_neo4j_driver, id="nonexistent")

        assert result is None

    @pytest.mark.asyncio
    async def test_get_multiple_ids(self, valid_node_class, mock_neo4j_driver):
        """Test get() with multiple ids parameter."""
        mock_records = [
            {"node": {"id": "id1", "name": "Node 1", "age": 10}},
            {"node": {"id": "id2", "name": "Node 2", "age": 20}},
        ]

        async def mock_execute_read(driver, query, params, processor):
            return mock_records

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            results = await valid_node_class.get(mock_neo4j_driver, ids=["id1", "id2"])

        assert len(results) == 2
        assert results[0].id == "id1"
        assert results[1].id == "id2"

    @pytest.mark.asyncio
    async def test_get_raises_error_both_params(self, valid_node_class, mock_neo4j_driver):
        """Test get() raises error when both id and ids provided."""
        with pytest.raises(ValueError, match="Provide either 'id' or 'ids', not both"):
            await valid_node_class.get(mock_neo4j_driver, id="test", ids=["test"])

    @pytest.mark.asyncio
    async def test_get_raises_error_no_params(self, valid_node_class, mock_neo4j_driver):
        """Test get() raises error when neither id nor ids provided."""
        with pytest.raises(ValueError, match="Provide either 'id' or 'ids'"):
            await valid_node_class.get(mock_neo4j_driver)

    @pytest.mark.asyncio
    async def test_get_by_valid_filters(self, valid_node_class, mock_neo4j_driver):
        """Test get_by() with valid field filters."""
        mock_records = [{"node": {"id": "test-123", "name": "Test Node", "age": 42}}]

        async def mock_execute_read(driver, query, params, processor):
            return mock_records

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            results = await valid_node_class.get_by(mock_neo4j_driver, name="Test Node")

        assert len(results) == 1
        assert results[0].name == "Test Node"

    @pytest.mark.asyncio
    async def test_get_by_invalid_filter_key(self, valid_node_class, mock_neo4j_driver):
        """Test get_by() raises error for invalid filter keys."""
        with pytest.raises(ValueError, match="Invalid filter keys"):
            await valid_node_class.get_by(mock_neo4j_driver, invalid_field="value")

    @pytest.mark.asyncio
    async def test_get_by_no_filters(self, valid_node_class, mock_neo4j_driver):
        """Test get_by() raises error when no filters provided."""
        with pytest.raises(ValueError, match="At least one filter or range must be provided"):
            await valid_node_class.get_by(mock_neo4j_driver)

    @pytest.mark.asyncio
    async def test_get_all_with_pagination(self, valid_node_class, mock_neo4j_driver):
        """Test get_all() with pagination."""
        mock_records = [
            {"node": {"id": "id1", "name": "Node 1"}},
            {"node": {"id": "id2", "name": "Node 2"}},
        ]

        async def mock_execute_read(driver, query, params, processor):
            return mock_records

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            results = await valid_node_class.get_all(mock_neo4j_driver, limit=10, offset=0)

        assert len(results) == 2

    @pytest.mark.asyncio
    async def test_count(self, valid_node_class, mock_neo4j_driver):
        """Test count() method."""
        mock_data = {"count": 42}

        async def mock_execute_read(driver, query, params, processor):
            return mock_data

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            count = await valid_node_class.count(mock_neo4j_driver)

        assert count == 42

    @pytest.mark.asyncio
    async def test_count_returns_zero_when_no_data(self, valid_node_class, mock_neo4j_driver):
        """Test count() returns 0 when no data returned."""

        async def mock_execute_read(driver, query, params, processor):
            return None

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            count = await valid_node_class.count(mock_neo4j_driver)

        assert count == 0


# ============================================================================
# Test get_by() with Ranges Parameter
# ============================================================================


class TestGetByWithRanges:
    """Tests for get_by() method with ranges parameter."""

    @pytest.fixture
    def numeric_node_class(self):
        """Create a node class with numeric fields for range testing."""

        class NumericNode(BaseNode):
            id: Annotated[str, LabelProperty(index=True)] = Field(..., description="ID")
            name: str = Field(..., description="Name")
            age: int | None = Field(None, description="Age")
            weight: float | None = Field(None, description="Weight")

        return NumericNode

    @pytest.mark.asyncio
    async def test_get_by_with_ranges_only(self, numeric_node_class, mock_neo4j_driver):
        """Test get_by() with ranges parameter only."""
        mock_records = [
            {"node": {"id": "id1", "name": "Node 1", "age": 25, "weight": 70.5}},
            {"node": {"id": "id2", "name": "Node 2", "age": 30, "weight": 75.0}},
        ]

        captured_query = None
        captured_params = None

        async def mock_execute_read(driver, query, params, processor):
            nonlocal captured_query, captured_params
            captured_query = query
            captured_params = params
            return mock_records

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            results = await numeric_node_class.get_by(
                mock_neo4j_driver,
                ranges={"age": (20, 35)},
            )

        assert len(results) == 2
        assert captured_query is not None
        assert "n.`age` >= $age_min" in captured_query
        assert "n.`age` <= $age_max" in captured_query
        assert captured_params["age_min"] == 20
        assert captured_params["age_max"] == 35

    @pytest.mark.asyncio
    async def test_get_by_with_ranges_and_exact_filters(self, numeric_node_class, mock_neo4j_driver):
        """Test get_by() with both ranges and exact filters."""
        mock_records = [{"node": {"id": "id1", "name": "Test", "age": 25, "weight": 70.5}}]

        captured_query = None
        captured_params = None

        async def mock_execute_read(driver, query, params, processor):
            nonlocal captured_query, captured_params
            captured_query = query
            captured_params = params
            return mock_records

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            results = await numeric_node_class.get_by(
                mock_neo4j_driver,
                name="Test",
                ranges={"age": (20, 35)},
            )

        assert len(results) == 1
        assert captured_query is not None
        # Verify both exact match and range conditions
        assert "n.`name` = $name" in captured_query
        assert "n.`age` >= $age_min" in captured_query
        assert "n.`age` <= $age_max" in captured_query
        assert captured_params["name"] == "Test"
        assert captured_params["age_min"] == 20
        assert captured_params["age_max"] == 35

    @pytest.mark.asyncio
    async def test_get_by_with_open_range_min_only(self, numeric_node_class, mock_neo4j_driver):
        """Test get_by() with range having only min value (no upper bound)."""
        mock_records = [{"node": {"id": "id1", "name": "Node 1", "age": 50, "weight": 80.0}}]

        captured_query = None
        captured_params = None

        async def mock_execute_read(driver, query, params, processor):
            nonlocal captured_query, captured_params
            captured_query = query
            captured_params = params
            return mock_records

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            results = await numeric_node_class.get_by(
                mock_neo4j_driver,
                ranges={"age": (30, None)},
            )

        assert len(results) == 1
        assert captured_query is not None
        assert "n.`age` >= $age_min" in captured_query
        assert "age_max" not in captured_query
        assert captured_params["age_min"] == 30
        assert "age_max" not in captured_params

    @pytest.mark.asyncio
    async def test_get_by_with_open_range_max_only(self, numeric_node_class, mock_neo4j_driver):
        """Test get_by() with range having only max value (no lower bound)."""
        mock_records = [{"node": {"id": "id1", "name": "Node 1", "age": 20, "weight": 60.0}}]

        captured_query = None
        captured_params = None

        async def mock_execute_read(driver, query, params, processor):
            nonlocal captured_query, captured_params
            captured_query = query
            captured_params = params
            return mock_records

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            results = await numeric_node_class.get_by(
                mock_neo4j_driver,
                ranges={"age": (None, 25)},
            )

        assert len(results) == 1
        assert captured_query is not None
        assert "age_min" not in captured_query
        assert "n.`age` <= $age_max" in captured_query
        assert "age_min" not in captured_params
        assert captured_params["age_max"] == 25

    @pytest.mark.asyncio
    async def test_get_by_with_multiple_ranges(self, numeric_node_class, mock_neo4j_driver):
        """Test get_by() with multiple range filters."""
        mock_records = [{"node": {"id": "id1", "name": "Node 1", "age": 30, "weight": 75.0}}]

        captured_query = None
        captured_params = None

        async def mock_execute_read(driver, query, params, processor):
            nonlocal captured_query, captured_params
            captured_query = query
            captured_params = params
            return mock_records

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            results = await numeric_node_class.get_by(
                mock_neo4j_driver,
                ranges={
                    "age": (25, 35),
                    "weight": (70.0, 80.0),
                },
            )

        assert len(results) == 1
        assert captured_query is not None
        # Verify both range conditions
        assert "n.`age` >= $age_min" in captured_query
        assert "n.`age` <= $age_max" in captured_query
        assert "n.`weight` >= $weight_min" in captured_query
        assert "n.`weight` <= $weight_max" in captured_query
        assert captured_params["age_min"] == 25
        assert captured_params["age_max"] == 35
        assert captured_params["weight_min"] == 70.0
        assert captured_params["weight_max"] == 80.0

    @pytest.mark.asyncio
    async def test_get_by_invalid_range_key(self, numeric_node_class, mock_neo4j_driver):
        """Test get_by() raises error for invalid range field names."""
        with pytest.raises(ValueError, match="Invalid range keys"):
            await numeric_node_class.get_by(
                mock_neo4j_driver,
                ranges={"invalid_field": (10, 20)},
            )

    @pytest.mark.asyncio
    async def test_get_by_empty_ranges_and_filters_raises_error(
        self, numeric_node_class, mock_neo4j_driver
    ):
        """Test get_by() raises error when both ranges and filters are empty."""
        with pytest.raises(ValueError, match="At least one filter or range must be provided"):
            await numeric_node_class.get_by(mock_neo4j_driver, ranges={})

    @pytest.mark.asyncio
    async def test_get_by_none_ranges_with_filters(self, numeric_node_class, mock_neo4j_driver):
        """Test get_by() works with ranges=None and exact filters provided."""
        mock_records = [{"node": {"id": "id1", "name": "Test", "age": 30, "weight": 70.0}}]

        async def mock_execute_read(driver, query, params, processor):
            return mock_records

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            results = await numeric_node_class.get_by(
                mock_neo4j_driver,
                name="Test",
                ranges=None,
            )

        assert len(results) == 1
        assert results[0].name == "Test"

    @pytest.mark.asyncio
    async def test_get_by_backward_compatibility(self, numeric_node_class, mock_neo4j_driver):
        """Test get_by() maintains backward compatibility without ranges parameter."""
        mock_records = [{"node": {"id": "id1", "name": "Test", "age": 30, "weight": 70.0}}]

        async def mock_execute_read(driver, query, params, processor):
            return mock_records

        with patch("pyeed.ingest.model.pyeedbase.execute_read", side_effect=mock_execute_read):
            # Call without ranges parameter (backward compatible)
            results = await numeric_node_class.get_by(mock_neo4j_driver, name="Test")

        assert len(results) == 1
        assert results[0].name == "Test"
