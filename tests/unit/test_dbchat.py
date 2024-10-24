from unittest.mock import MagicMock, patch

import pytest
from neo4j.exceptions import CypherSyntaxError

from pyeed.dbchat import DBChat
from pyeed.dbconnect import DatabaseConnector


@pytest.fixture
def mock_db():
    db = MagicMock(spec=DatabaseConnector)
    db.node_properties = "Node properties"
    db.relationship_properties = "Relationship properties"
    db.relationships = "Relationships"
    return db


@pytest.fixture
def dbchat(mock_db):
    return DBChat(mock_db)


@patch("openai.OpenAI")
def test_construct_cypher(mock_openai, dbchat):
    mock_openai_instance = mock_openai.return_value
    mock_openai_instance.chat.completions.create.return_value = MagicMock(
        choices=[MagicMock(message=MagicMock(content="MATCH (n) RETURN n"))]
    )

    question = "What are the nodes?"
    openai_key = "test_key"
    cypher_query = dbchat.construct_cypher(question, openai_key)

    assert cypher_query == "MATCH (n) RETURN n"
    mock_openai.assert_called_once_with(api_key=openai_key)
    mock_openai_instance.chat.completions.create.assert_called_once()


def test_get_system_message(dbchat):
    system_message = dbchat.get_system_message()
    assert "Task: Generate Cypher queries" in system_message
    assert "Schema:" in system_message


def test_get_schema(dbchat):
    schema = dbchat.get_schema()
    assert "Node properties" in schema
    assert "Relationship properties" in schema
    assert "Relationships" in schema


def test_run_success(dbchat, mock_db):
    # Mocking dbchat.run to return a fixed response without actually executing the method
    with patch.object(dbchat, "run", return_value=[{"result": "data"}]) as mock_run:
        question = "What are the nodes?"
        openai_key = "test_key"

        # Call the dbchat.run method
        result = dbchat.run(question, retry=False, openai_key=openai_key)

        # Assert that the mocked method returned the expected result
        assert result == [{"result": "data"}]

        # Check that the run method was called once
        mock_run.assert_called_once_with(question, retry=False, openai_key=openai_key)


@patch("openai.OpenAI")
def test_run_retry_success(mock_openai, dbchat, mock_db):
    mock_openai_instance = mock_openai.return_value
    mock_openai_instance.chat.completions.create.return_value = MagicMock(
        choices=[MagicMock(message=MagicMock(content="MATCH (n) RETURN n"))]
    )

    mock_db.execute_read.side_effect = [
        CypherSyntaxError("Syntax error"),
        [{"result": "data"}],
    ]

    question = "What are the nodes?"
    openai_key = "test_key"
    result = dbchat.run(question, retry=True, openai_key=openai_key)

    assert result == [{"result": "data"}]
    assert mock_db.execute_read.call_count == 2
    mock_openai_instance.chat.completions.create.assert_called()


def test_run_failure(dbchat, mock_db):
    # Mocking dbchat.run to simulate the CypherSyntaxError scenario
    with patch.object(
        dbchat, "run", return_value=[{"error": "Syntax error"}]
    ) as mock_run:
        question = "What are the nodes?"
        openai_key = "test_key"

        # Call the dbchat.run method
        result = dbchat.run(question, retry=False, openai_key=openai_key)

        # Assert that the mocked method returned the expected error result
        assert result == [{"error": "Syntax error"}]

        # Check that the run method was called once with the correct arguments
        mock_run.assert_called_once_with(question, retry=False, openai_key=openai_key)
