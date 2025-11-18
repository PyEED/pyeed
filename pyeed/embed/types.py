"""Type definitions for embeddings."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal, Protocol

import numpy as np
import torch
from numpy.typing import NDArray

if TYPE_CHECKING:
    from ..ingest.core.pipeline import PipelineRecord
    from ..ingest.model.pyeedbase import PyeedBase

__all__ = [
    "NP_DTYPE_MAP",
    "TF_DTYPE_MAP",
    "EmbeddableData",
    "EmbeddingBatch",
    "EmbeddingRecord",
    "ModelDType",
    "ReturnDType",
    "Vector",
    "distribute_embeddings_to_records",
    "flatten_embedding_record",
    "records_to_embedding_inputs",
]

# Type literals for model and return dtypes
type ModelDType = Literal["float16", "bfloat16", "float32"]
type ReturnDType = Literal["float16", "bfloat16", "float32"]


type Vector = NDArray[np.float32] | NDArray[np.float16]


# Dtype mappings
TF_DTYPE_MAP: dict[ModelDType, torch.dtype] = {
    "float16": torch.float16,
    "float32": torch.float32,
}

NP_DTYPE_MAP: dict[ReturnDType, type[np.generic]] = {
    "float16": np.float16,
    "float32": np.float32,
}


class EmbeddableData(Protocol):
    """Protocol for data objects that can be embedded.

    Any object with id and sequence attributes can be converted to embedder inputs.
    """

    id: str
    sequence: str


@dataclass(slots=True)
class EmbeddingBatch:
    protein_ids: list[str] = field(default_factory=list)
    sequences: list[str] = field(default_factory=list)
    embeddings: dict[str, list[Vector]] = field(default_factory=dict)


@dataclass(slots=True)
class EmbeddingRecord:
    protein_id: str
    sequence: str
    embedding: dict[str, Vector]


def flatten_embedding_record(
    record: EmbeddingRecord,
) -> dict[str, str | Vector]:
    """Flatten an embedding record into a dictionary."""

    d: dict[str, str | Vector] = {"protein_id": record.protein_id, "sequence": record.sequence}
    d.update(record.embedding)
    return d


def records_to_embedding_inputs[T: PyeedBase](
    records: list[PipelineRecord[T]],
) -> tuple[list[str], list[str]]:
    """Extract sequences and protein IDs from pipeline records for embedder input.

    Args:
        records: Pipeline records with embeddable data (must have id and sequence attributes)

    Returns:
        Tuple of (sequences, protein_ids) ready for embedder.embed_batch()

    Raises:
        AttributeError: If record.data doesn't have id or sequence attributes

    Example:
        >>> records = [PipelineRecord(data=protein) for protein in proteins]
        >>> sequences, protein_ids = records_to_embedding_inputs(records)
        >>> async for batch in embedder.embed_batch(sequences, protein_ids, batch_size=32):
        ...     updated = distribute_embeddings_to_records(batch, records)
    """
    sequences = [record.data.sequence for record in records]
    protein_ids = [record.data.id for record in records]
    return sequences, protein_ids


def distribute_embeddings_to_records[T: PyeedBase](
    embedding_batch: EmbeddingBatch,
    records: list[PipelineRecord[T]],
) -> list[PipelineRecord[T]]:
    """Distribute embeddings from batch back to original pipeline records.

    Maps embeddings to records by matching protein_ids. Updates record.embeddings
    dict in-place with all pooling methods from the batch.

    Args:
        embedding_batch: Batch of embeddings from embedder
        records: Original pipeline records to update

    Returns:
        Updated records with embeddings populated (same objects, modified in-place)

    Raises:
        ValueError: If protein_id from batch is not found in records

    Example:
        >>> sequences, protein_ids = records_to_embedding_inputs(records)
        >>> async for batch in embedder.embed_batch(sequences, protein_ids, batch_size=32):
        ...     updated = distribute_embeddings_to_records(batch, records)
        ...     for record in updated:
        ...         await output_queue.put(record)
    """
    # Build lookup: protein_id -> record
    record_map: dict[str, PipelineRecord[T]] = {record.data.id: record for record in records}

    # Distribute embeddings by matching protein_ids
    for idx, protein_id in enumerate(embedding_batch.protein_ids):
        if protein_id not in record_map:
            raise ValueError(
                f"Protein ID '{protein_id}' from embedding batch not found in records. "
                f"Ensure batch was created from these records."
            )

        record = record_map[protein_id]

        # Add all pooling methods to this record
        for pooling_method, embedding_list in embedding_batch.embeddings.items():
            record.embeddings[pooling_method] = embedding_list[idx]

    return records
