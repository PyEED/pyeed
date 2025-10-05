"""
DEPRECATED: Old embedding implementation.

This file contains the previous embedding system.
Use the new async API from embeddings/__init__.py instead.
"""

from typing import Any, Literal, cast

import numpy as np
import torch
from numpy.typing import NDArray
from torch.nn import DataParallel, Module
from transformers import EsmModel, EsmTokenizer, T5Model, T5Tokenizer

from pyeed.dbconnect import DatabaseConnector

from .base_deprecated import BaseEmbeddingModel, ModelType, normalize_embedding
from .database_deprecated import update_protein_embeddings_in_db
from .factory_deprecated import ModelFactory
from .models_deprecated import (
    ESM2EmbeddingModel,
    ESM3EmbeddingModel,
    ESMCEmbeddingModel,
    ProtT5EmbeddingModel,
)
from .processor_deprecated import EmbeddingProcessor, get_processor
from .utils import (
    determine_model_type,
    free_memory,
    get_hf_token,
    preprocess_sequence_for_prott5,
)


def load_model_and_tokenizer(
    model_name: str,
    device: torch.device | None = None,
) -> tuple[Any, Any | None, torch.device]:
    """DEPRECATED: Load model and tokenizer."""
    if device is None:
        device = torch.device("cuda:0")
    return ModelFactory.load_model_and_tokenizer(model_name, device)


def process_batches_on_gpu(
    data: list[tuple[str, str]],
    batch_size: int,
    model: Any | DataParallel[Module],
    tokenizer: Any | None,
    db: DatabaseConnector,
    device: torch.device,
) -> None:
    """DEPRECATED: Process batches on GPU."""
    processor = get_processor()
    processor.process_batches_on_gpu(data, batch_size, model, tokenizer, db, device)


def get_batch_embeddings(
    batch_sequences: list[str],
    model: Any | DataParallel[Module],
    tokenizer_or_alphabet: Any | None,
    device: torch.device,
    pool_embeddings: bool = True,
) -> list[NDArray[np.float64]]:
    """DEPRECATED: Get batch embeddings."""
    processor = get_processor()
    return processor.get_batch_embeddings_unified(
        batch_sequences, model, tokenizer_or_alphabet, device, pool_embeddings
    )


def calculate_single_sequence_embedding_last_hidden_state(
    sequence: str,
    device: torch.device | None = None,
    model_name: str = "facebook/esm2_t33_650M_UR50D",
) -> NDArray[np.float64]:
    """DEPRECATED: Calculate single sequence embedding using last hidden state."""
    if device is None:
        device = torch.device("cuda:0")
    processor = get_processor()
    return processor.calculate_single_sequence_embedding_last_hidden_state(
        sequence, device, model_name
    )


__all__ = [
    "load_model_and_tokenizer",
    "process_batches_on_gpu",
    "get_batch_embeddings",
    "calculate_single_sequence_embedding_last_hidden_state",
]


