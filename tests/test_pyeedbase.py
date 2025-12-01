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
# Test Custom Field Validation
# ============================================================================


class TestCustomFieldValidation:
    """Tests for custom field validation."""

    def test_valid_custom_fields(self, valid_node_class):
        """Test that valid custom fields are accepted."""
        node = valid_node_class(
            id="test",
            name="Node",
            custom={"extra": "value", "count": 42, "tags": ["a", "b"], "flag": True},
        )
        assert node.custom["extra"] == "value"
        assert node.custom["count"] == 42
        assert node.custom["tags"] == ["a", "b"]
        assert node.custom["flag"] is True

    def test_dict_value_raises_error(self, valid_node_class):
        """Test that dict values in custom fields raise ValueError."""
        with pytest.raises(ValueError, match="Nested dictionaries are not allowed"):
            valid_node_class(id="test", name="Node", custom={"nested": {"key": "value"}})

    def test_conflicting_key_raises_error(self, valid_node_class):
        """Test that keys conflicting with field names raise ValueError."""
        with pytest.raises(ValueError, match="cannot conflict with existing attributes"):
            valid_node_class(id="test", name="Node", custom={"id": "conflict"})

        with pytest.raises(ValueError, match="cannot conflict with existing attributes"):
            valid_node_class(id="test", name="Node", custom={"name": "conflict"})

    def test_custom_key_raises_error(self, valid_node_class):
        """Test that 'custom' as a key raises ValueError."""
        with pytest.raises(ValueError, match="cannot conflict with existing attributes"):
            valid_node_class(id="test", name="Node", custom={"custom": "value"})

    def test_invalid_key_format_raises_error(self, valid_node_class):
        """Test that invalid Python variable names raise ValueError."""
        with pytest.raises(ValueError, match="Invalid custom field keys"):
            valid_node_class(id="test", name="Node", custom={"123invalid": "value"})

        with pytest.raises(ValueError, match="Invalid custom field keys"):
            valid_node_class(id="test", name="Node", custom={"my-key": "value"})

        with pytest.raises(ValueError, match="Invalid custom field keys"):
            valid_node_class(id="test", name="Node", custom={"my key": "value"})

    def test_empty_custom_allowed(self, valid_node_class):
        """Test that empty custom dict is allowed."""
        node = valid_node_class(id="test", name="Node", custom={})
        assert node.custom == {}

    def test_none_values_in_custom_allowed(self, valid_node_class):
        """Test that None values in custom fields are allowed."""
        node = valid_node_class(id="test", name="Node", custom={"optional": None})
        assert node.custom["optional"] is None


# ============================================================================
# Test to_dict Method
# ============================================================================


class TestToDictMethod:
    """Tests for to_dict() method."""

    def test_basic_to_dict(self, sample_node):
        """Test basic to_dict conversion."""
        d = sample_node.to_dict()
        assert d["id"] == "test-123"
        assert d["name"] == "Test Node"
        assert d["age"] == 42

    def test_custom_fields_flattened(self, valid_node_class):
        """Test that custom fields are flattened into the dict."""
        node = valid_node_class(id="test", name="Node", custom={"extra": "value", "count": 42})
        d = node.to_dict()
        assert d["extra"] == "value"
        assert d["count"] == 42
        assert "custom" not in d

    def test_none_values_excluded(self, valid_node_class):
        """Test that None values are excluded from dict."""
        node = valid_node_class(id="test", name="Node", age=None)
        d = node.to_dict()
        assert "age" not in d

    def test_only_neo4j_compatible_values(self, valid_node_class):
        """Test that only Neo4j-compatible values are included."""
        node = valid_node_class(id="test", name="Node")
        # Add a complex object that shouldn't be included
        node.custom = {"valid_str": "text", "valid_list": [1, 2, 3]}
        d = node.to_dict()
        assert "valid_str" in d
        assert "valid_list" in d

    def test_list_of_primitives_included(self, valid_node_class):
        """Test that lists of primitives are included."""
        node = valid_node_class(
            id="test", name="Node", custom={"tags": ["a", "b"], "numbers": [1, 2, 3]}
        )
        d = node.to_dict()
        assert d["tags"] == ["a", "b"]
        assert d["numbers"] == [1, 2, 3]


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
    """Tests for bulk_upsert() method with mocked Neo4j."""

    @pytest.mark.asyncio
    async def test_bulk_upsert_basic(self, valid_node_class, mock_neo4j_driver):
        """Test basic bulk upsert operation."""
        nodes = [
            valid_node_class(id="id1", name="Node 1"),
            valid_node_class(id="id2", name="Node 2"),
        ]

        await valid_node_class.bulk_upsert(mock_neo4j_driver, nodes)

        # Verify driver.session was called
        mock_neo4j_driver.session.assert_called_once()

    @pytest.mark.asyncio
    async def test_bulk_upsert_with_session_kwargs(self, valid_node_class, mock_neo4j_driver):
        """Test bulk upsert with session kwargs."""
        nodes = [valid_node_class(id="id1", name="Node 1")]
        session_kwargs = {"database": "test_db"}

        await valid_node_class.bulk_upsert(mock_neo4j_driver, nodes, session_kwargs=session_kwargs)

        # Verify session was called with kwargs
        mock_neo4j_driver.session.assert_called_once_with(**session_kwargs)

    @pytest.mark.asyncio
    async def test_bulk_upsert_empty_list(self, valid_node_class, mock_neo4j_driver):
        """Test bulk upsert with empty list does nothing."""
        await valid_node_class.bulk_upsert(mock_neo4j_driver, [])

        # Should not call session
        mock_neo4j_driver.session.assert_not_called()

    @pytest.mark.asyncio
    async def test_bulk_upsert_generates_correct_cypher(self, valid_node_class, mock_neo4j_driver):
        """Test that bulk upsert generates correct Cypher query."""
        nodes = [valid_node_class(id="id1", name="Node 1")]

        # Mock the transaction execution to capture the cypher
        captured_cypher = None
        captured_rows = None

        async def mock_execute_write(tx_func, rows):
            nonlocal captured_cypher, captured_rows
            # Create a mock transaction
            mock_tx = AsyncMock()

            async def capture_run(cypher, **kwargs):
                nonlocal captured_cypher, captured_rows
                captured_cypher = cypher
                captured_rows = kwargs.get("rows")

            mock_tx.run = capture_run
            await tx_func(mock_tx, rows)

        mock_neo4j_driver.session.return_value.__aenter__.return_value.execute_write = (
            mock_execute_write
        )

        await valid_node_class.bulk_upsert(mock_neo4j_driver, nodes)

        # Verify Cypher contains expected elements
        assert captured_cypher is not None
        assert "UNWIND $rows AS row" in captured_cypher
        assert "MERGE" in captured_cypher
        assert "TestNode" in captured_cypher
        assert "id: row.id" in captured_cypher
        assert "SET n += row" in captured_cypher

        # Verify rows data
        assert captured_rows is not None
        assert len(captured_rows) == 1
        assert captured_rows[0]["id"] == "id1"
        assert captured_rows[0]["name"] == "Node 1"

    @pytest.mark.asyncio
    async def test_bulk_upsert_with_custom_fields(self, valid_node_class, mock_neo4j_driver):
        """Test bulk upsert includes custom fields."""
        nodes = [valid_node_class(id="id1", name="Node 1", custom={"extra": "value", "count": 42})]

        captured_rows = None

        async def mock_execute_write(tx_func, rows):
            nonlocal captured_rows
            captured_rows = rows
            mock_tx = AsyncMock()
            mock_tx.run = AsyncMock()
            await tx_func(mock_tx, rows)

        mock_neo4j_driver.session.return_value.__aenter__.return_value.execute_write = (
            mock_execute_write
        )

        await valid_node_class.bulk_upsert(mock_neo4j_driver, nodes)

        # Verify custom fields are flattened into rows
        assert captured_rows[0]["extra"] == "value"
        assert captured_rows[0]["count"] == 42


class TestUpsert:
    """Tests for single node upsert() method."""

    @pytest.mark.asyncio
    async def test_upsert_calls_bulk_upsert(self, sample_node, mock_neo4j_driver):
        """Test that upsert calls bulk_upsert with single node."""
        with patch.object(type(sample_node), "bulk_upsert", new_callable=AsyncMock) as mock_bulk:
            await sample_node.upsert(mock_neo4j_driver)

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

        with patch.object(type(sample_node), "bulk_upsert", new_callable=AsyncMock) as mock_bulk:
            await sample_node.upsert(mock_neo4j_driver, session_kwargs=session_kwargs)

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
        # Create node with custom fields
        node = valid_node_class(
            id="test-123",
            name="Integration Test",
            age=99,
            custom={"tags": ["tag1", "tag2"], "score": 95.5},
        )

        # Verify node creation
        assert node.id == "test-123"
        assert node.name == "Integration Test"
        assert node.age == 99

        # Verify get_index_field
        assert node.get_index_field() == "id"

        # Verify to_dict includes everything
        d = node.to_dict()
        assert d["id"] == "test-123"
        assert d["name"] == "Integration Test"
        assert d["age"] == 99
        assert d["tags"] == ["tag1", "tag2"]
        assert d["score"] == 95.5
        assert "custom" not in d

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

        # Verify to_dict includes all fields
        d = child.to_dict()
        assert d["id"] == "child-1"
        assert d["parent_field"] == "parent"
        assert d["child_field"] == "child"
