"""Unit tests for embedding type conversion utilities."""

import numpy as np
import pytest

from pyeed.embed.types import (
    EmbeddingBatch,
    distribute_embeddings_to_records,
    records_to_embedding_inputs,
)
from pyeed.ingest.core.pipeline import PipelineRecord
from pyeed.ingest.model.protein import Protein


@pytest.fixture
def sample_proteins() -> list[Protein]:
    """Create sample protein objects for testing."""
    return [
        Protein(
            id="P12345",
            sequence="MKTAYIAKQRQISFVKSHFS",
            seq_length=20,
            name="Test Protein 1",
        ),
        Protein(
            id="Q67890",
            sequence="ARNDCEQGHILKMFPSTWYV",
            seq_length=20,
            name="Test Protein 2",
        ),
        Protein(
            id="R11111",
            sequence="MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDLLKRMKV",
            seq_length=49,
            name="Test Protein 3",
        ),
    ]


@pytest.fixture
def pipeline_records(sample_proteins: list[Protein]) -> list[PipelineRecord[Protein]]:
    """Create pipeline records from sample proteins."""
    return [PipelineRecord(data=protein) for protein in sample_proteins]


@pytest.fixture
def embedding_batch(sample_proteins: list[Protein]) -> EmbeddingBatch:
    """Create a mock embedding batch with multiple pooling methods."""
    batch = EmbeddingBatch()
    batch.protein_ids = [p.id for p in sample_proteins]
    batch.sequences = [p.sequence for p in sample_proteins]

    # Create mock embeddings for two pooling methods
    embedding_dim = 128
    batch.embeddings = {
        "mean_pool": [np.random.randn(embedding_dim).astype(np.float32) for _ in sample_proteins],
        "cls_token": [np.random.randn(embedding_dim).astype(np.float32) for _ in sample_proteins],
    }

    return batch


class TestRecordsToEmbeddingInputs:
    """Tests for records_to_embedding_inputs conversion function."""

    def test_basic_conversion(self, pipeline_records: list[PipelineRecord[Protein]]) -> None:
        """Test basic conversion of records to embedder inputs."""
        sequences, protein_ids = records_to_embedding_inputs(pipeline_records)

        assert len(sequences) == len(pipeline_records)
        assert len(protein_ids) == len(pipeline_records)

        for i, record in enumerate(pipeline_records):
            assert sequences[i] == record.data.sequence
            assert protein_ids[i] == record.data.id

    def test_empty_records(self) -> None:
        """Test conversion with empty record list."""
        sequences, protein_ids = records_to_embedding_inputs([])

        assert sequences == []
        assert protein_ids == []

    def test_single_record(self, sample_proteins: list[Protein]) -> None:
        """Test conversion with single record."""
        record = PipelineRecord(data=sample_proteins[0])
        sequences, protein_ids = records_to_embedding_inputs([record])

        assert len(sequences) == 1
        assert len(protein_ids) == 1
        assert sequences[0] == sample_proteins[0].sequence
        assert protein_ids[0] == sample_proteins[0].id

    def test_preserves_order(self, pipeline_records: list[PipelineRecord[Protein]]) -> None:
        """Test that order is preserved in conversion."""
        sequences, protein_ids = records_to_embedding_inputs(pipeline_records)

        for i, record in enumerate(pipeline_records):
            assert sequences[i] == record.data.sequence
            assert protein_ids[i] == record.data.id


class TestDistributeEmbeddingsToRecords:
    """Tests for distribute_embeddings_to_records conversion function."""

    def test_basic_distribution(
        self,
        pipeline_records: list[PipelineRecord[Protein]],
        embedding_batch: EmbeddingBatch,
    ) -> None:
        """Test basic distribution of embeddings to records."""
        updated_records = distribute_embeddings_to_records(embedding_batch, pipeline_records)

        # Should return the same record objects
        assert updated_records is pipeline_records

        # Check each record has embeddings
        expected_pooling_methods = 2  # mean_pool and cls_token
        embedding_dim = 128
        for record in updated_records:
            assert len(record.embeddings) == expected_pooling_methods
            assert "mean_pool" in record.embeddings
            assert "cls_token" in record.embeddings

            # Verify embedding shapes
            assert record.embeddings["mean_pool"].shape == (embedding_dim,)
            assert record.embeddings["cls_token"].shape == (embedding_dim,)

    def test_protein_id_matching(
        self,
        pipeline_records: list[PipelineRecord[Protein]],
        embedding_batch: EmbeddingBatch,
    ) -> None:
        """Test that embeddings are correctly matched by protein_id."""
        # Shuffle the records order
        shuffled_records = [pipeline_records[2], pipeline_records[0], pipeline_records[1]]

        updated_records = distribute_embeddings_to_records(embedding_batch, shuffled_records)

        # Each record should still get the correct embedding
        for record in updated_records:
            # Find the index in the original batch
            batch_idx = embedding_batch.protein_ids.index(record.data.id)

            # Verify the embeddings match
            for pooling_method, embeddings in embedding_batch.embeddings.items():
                np.testing.assert_array_equal(
                    record.embeddings[pooling_method],
                    embeddings[batch_idx],
                )

    def test_missing_protein_id_raises_error(
        self,
        pipeline_records: list[PipelineRecord[Protein]],
        embedding_batch: EmbeddingBatch,
    ) -> None:
        """Test that missing protein_id raises ValueError."""
        # Add a protein_id to batch that doesn't exist in records
        embedding_batch.protein_ids.append("MISSING_ID")
        embedding_batch.embeddings["mean_pool"].append(np.zeros(128, dtype=np.float32))
        embedding_batch.embeddings["cls_token"].append(np.zeros(128, dtype=np.float32))

        with pytest.raises(ValueError, match="Protein ID 'MISSING_ID' from embedding batch"):
            distribute_embeddings_to_records(embedding_batch, pipeline_records)

    def test_multiple_pooling_methods(
        self,
        pipeline_records: list[PipelineRecord[Protein]],
        sample_proteins: list[Protein],
    ) -> None:
        """Test distribution with multiple pooling methods."""
        # Create batch with 3 different pooling methods
        batch = EmbeddingBatch()
        batch.protein_ids = [p.id for p in sample_proteins]
        batch.sequences = [p.sequence for p in sample_proteins]
        batch.embeddings = {
            "mean_pool": [np.ones(64, dtype=np.float32) * i for i in range(3)],
            "max_pool": [np.ones(64, dtype=np.float32) * (i + 10) for i in range(3)],
            "cls_token": [np.ones(64, dtype=np.float32) * (i + 20) for i in range(3)],
        }

        updated_records = distribute_embeddings_to_records(batch, pipeline_records)

        expected_pooling_methods = 3
        for i, record in enumerate(updated_records):
            assert len(record.embeddings) == expected_pooling_methods
            np.testing.assert_array_equal(record.embeddings["mean_pool"], np.ones(64) * i)
            np.testing.assert_array_equal(record.embeddings["max_pool"], np.ones(64) * (i + 10))
            np.testing.assert_array_equal(record.embeddings["cls_token"], np.ones(64) * (i + 20))

    def test_empty_batch(self) -> None:
        """Test distribution with empty batch."""
        empty_batch = EmbeddingBatch()
        empty_records: list[PipelineRecord[Protein]] = []

        updated_records = distribute_embeddings_to_records(empty_batch, empty_records)
        assert updated_records == []

    def test_modifies_in_place(
        self,
        pipeline_records: list[PipelineRecord[Protein]],
        embedding_batch: EmbeddingBatch,
    ) -> None:
        """Test that records are modified in-place."""
        original_ids = [id(record) for record in pipeline_records]

        updated_records = distribute_embeddings_to_records(embedding_batch, pipeline_records)

        # Verify same object references
        updated_ids = [id(record) for record in updated_records]
        assert original_ids == updated_ids

        # Verify embeddings were added to original objects
        for record in pipeline_records:
            assert len(record.embeddings) > 0


class TestRoundTrip:
    """Integration tests for round-trip conversion."""

    def test_full_pipeline_flow(
        self,
        pipeline_records: list[PipelineRecord[Protein]],
        embedding_batch: EmbeddingBatch,
    ) -> None:
        """Test complete flow: records -> inputs -> embeddings -> back to records."""
        # Step 1: Convert records to embedder inputs
        sequences, protein_ids = records_to_embedding_inputs(pipeline_records)

        # Step 2: Verify embedder inputs match what we'd expect
        assert sequences == [r.data.sequence for r in pipeline_records]
        assert protein_ids == [r.data.id for r in pipeline_records]

        # Step 3: Distribute embeddings back to records
        updated_records = distribute_embeddings_to_records(embedding_batch, pipeline_records)

        # Step 4: Verify all records have embeddings
        for record in updated_records:
            assert record.embeddings
            assert "mean_pool" in record.embeddings
            assert "cls_token" in record.embeddings
            assert isinstance(record.embeddings["mean_pool"], np.ndarray)
            assert isinstance(record.embeddings["cls_token"], np.ndarray)

    def test_multiple_batches(
        self,
        sample_proteins: list[Protein],
    ) -> None:
        """Test processing multiple batches and accumulating embeddings."""
        # Create two separate batches
        records_batch1 = [PipelineRecord(data=sample_proteins[0])]
        records_batch2 = [
            PipelineRecord(data=sample_proteins[1]),
            PipelineRecord(data=sample_proteins[2]),
        ]

        # Process batch 1
        sequences1, protein_ids1 = records_to_embedding_inputs(records_batch1)
        embed_batch1 = EmbeddingBatch(
            protein_ids=protein_ids1,
            sequences=sequences1,
            embeddings={"mean_pool": [np.ones(64, dtype=np.float32)]},
        )
        distribute_embeddings_to_records(embed_batch1, records_batch1)

        # Process batch 2
        sequences2, protein_ids2 = records_to_embedding_inputs(records_batch2)
        embed_batch2 = EmbeddingBatch(
            protein_ids=protein_ids2,
            sequences=sequences2,
            embeddings={
                "mean_pool": [np.ones(64, dtype=np.float32) * 2, np.ones(64, dtype=np.float32) * 3]
            },
        )
        distribute_embeddings_to_records(embed_batch2, records_batch2)

        # Verify each batch has correct embeddings
        expected_value_batch1 = 1.0
        expected_value_batch2_first = 2.0
        expected_value_batch2_second = 3.0
        assert records_batch1[0].embeddings["mean_pool"][0] == expected_value_batch1
        assert records_batch2[0].embeddings["mean_pool"][0] == expected_value_batch2_first
        assert records_batch2[1].embeddings["mean_pool"][0] == expected_value_batch2_second
