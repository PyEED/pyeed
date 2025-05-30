"""
Organized embedding module for protein language models.

This module provides both the new organized structure and backward compatibility
with the original embedding.py interface.
"""

from typing import Any, Tuple, Union, List, Optional, cast
import torch
from torch.nn import DataParallel, Module
from numpy.typing import NDArray
import numpy as np
from transformers import EsmModel, EsmTokenizer, T5Model, T5Tokenizer
from esm.models.esmc import ESMC
from esm.models.esm3 import ESM3

# New organized structure
from .base import BaseEmbeddingModel, ModelType, normalize_embedding
from .factory import ModelFactory
from .processor import EmbeddingProcessor, get_processor
from .utils import get_hf_token, preprocess_sequence_for_prott5, free_memory, determine_model_type
from .database import update_protein_embeddings_in_db
from .models import ESM2EmbeddingModel, ESMCEmbeddingModel, ESM3EmbeddingModel, ProtT5EmbeddingModel

from pyeed.dbconnect import DatabaseConnector

# Type aliases for better readability
TokenizerType = Union[EsmTokenizer, T5Tokenizer, None]
DeviceType = torch.device

# Re-export functions from processor
__all__ = [
    'load_model_and_tokenizer',
    'process_batches_on_gpu',
    'get_batch_embeddings',
    'calculate_single_sequence_embedding_last_hidden_state',
    'calculate_single_sequence_embedding_all_layers',
    'calculate_single_sequence_embedding_first_layer',
    'get_single_embedding_last_hidden_state',
    'get_single_embedding_all_layers',
    'get_single_embedding_first_layer',
]

# Function implementations
def load_model_and_tokenizer(
    model_name: str,
    device: Optional[DeviceType] = None,
) -> Tuple[ModelType, TokenizerType, DeviceType]:
    """Load model and tokenizer."""
    if device is None:
        device = torch.device("cuda:0")
    return cast(Tuple[ModelType, TokenizerType, DeviceType], ModelFactory.load_model_and_tokenizer(model_name, device))


def process_batches_on_gpu(
    data: List[Tuple[str, str]],
    batch_size: int,
    model: Union[EsmModel, ESMC, ESM3, T5Model, DataParallel[Module]],
    tokenizer: Union[EsmTokenizer, T5Tokenizer, None],
    db: DatabaseConnector,
    device: torch.device,
) -> None:
    """Process batches on GPU."""
    processor = get_processor()
    processor.process_batches_on_gpu(data, batch_size, model, tokenizer, db, device)


def get_batch_embeddings(
    batch_sequences: List[str],
    model: Union[EsmModel, ESMC, ESM3, T5Model, DataParallel[Module]],
    tokenizer_or_alphabet: Union[EsmTokenizer, T5Tokenizer, None],
    device: torch.device,
    pool_embeddings: bool = True,
) -> List[NDArray[np.float64]]:
    """Get batch embeddings."""
    processor = get_processor()
    return processor.get_batch_embeddings_unified(
        batch_sequences, model, tokenizer_or_alphabet, device, pool_embeddings
    )


def calculate_single_sequence_embedding_last_hidden_state(
    sequence: str,
    device: Optional[torch.device] = None,
    model_name: str = "facebook/esm2_t33_650M_UR50D",
) -> NDArray[np.float64]:
    """Calculate single sequence embedding using last hidden state."""
    if device is None:
        device = torch.device("cuda:0")
    processor = get_processor()
    return processor.calculate_single_sequence_embedding_last_hidden_state(
        sequence, device, model_name
    )


def calculate_single_sequence_embedding_all_layers(
    sequence: str,
    device: torch.device,
    model_name: str = "facebook/esm2_t33_650M_UR50D",
) -> NDArray[np.float64]:
    """Calculate single sequence embedding using all layers."""
    processor = get_processor()
    return processor.calculate_single_sequence_embedding_all_layers(
        sequence, device, model_name
    )


def calculate_single_sequence_embedding_first_layer(
    sequence: str,
    model_name: str = "facebook/esm2_t33_650M_UR50D",
    device: Optional[torch.device] = None,
) -> NDArray[np.float64]:
    """Calculate single sequence embedding using first layer."""
    if device is None:
        device = torch.device("cuda:0")
    processor = get_processor()
    return processor.calculate_single_sequence_embedding_first_layer(
        sequence, model_name, device
    )


def get_single_embedding_last_hidden_state(
    sequence: str,
    model: Union[EsmModel, ESMC, ESM3, T5Model, DataParallel[Module]],
    tokenizer: Union[EsmTokenizer, T5Tokenizer, None],
    device: torch.device,
) -> NDArray[np.float64]:
    """Get single embedding using last hidden state."""
    processor = get_processor()
    return processor.get_single_embedding_last_hidden_state(sequence, model, tokenizer, device)


def get_single_embedding_all_layers(
    sequence: str,
    model: Union[EsmModel, ESMC, ESM3, T5Model, DataParallel[Module]],
    tokenizer: Union[EsmTokenizer, T5Tokenizer, None],
    device: torch.device,
) -> NDArray[np.float64]:
    """Get single embedding using all layers."""
    processor = get_processor()
    return processor.get_single_embedding_all_layers(sequence, model, tokenizer, device)


def get_single_embedding_first_layer(
    sequence: str,
    model: Union[EsmModel, ESMC, ESM3, T5Model, DataParallel[Module]],
    tokenizer: Union[EsmTokenizer, T5Tokenizer, None],
    device: torch.device,
) -> NDArray[np.float64]:
    """Get single embedding using first layer."""
    processor = get_processor()
    return processor.get_single_embedding_first_layer(sequence, model, tokenizer, device)

# Public API
load_model_and_tokenizer = load_model_and_tokenizer
process_batches_on_gpu = process_batches_on_gpu
get_batch_embeddings = get_batch_embeddings
calculate_single_sequence_embedding_last_hidden_state = calculate_single_sequence_embedding_last_hidden_state
calculate_single_sequence_embedding_all_layers = calculate_single_sequence_embedding_all_layers
calculate_single_sequence_embedding_first_layer = calculate_single_sequence_embedding_first_layer
get_single_embedding_last_hidden_state = get_single_embedding_last_hidden_state
get_single_embedding_all_layers = get_single_embedding_all_layers
get_single_embedding_first_layer = get_single_embedding_first_layer

__all__ = [
    # Base classes and types
    'BaseEmbeddingModel',
    'ModelType',
    'normalize_embedding',
    
    # Factory and processor
    'ModelFactory',
    'EmbeddingProcessor',
    'get_processor',
    
    # Utilities
    'get_hf_token',
    'preprocess_sequence_for_prott5',
    'free_memory',
    'determine_model_type',
    
    # Database operations
    'update_protein_embeddings_in_db',
    
    # Model implementations
    'ESM2EmbeddingModel',
    'ESMCEmbeddingModel',
    'ESM3EmbeddingModel',
    'ProtT5EmbeddingModel',
    
    # Backward compatibility functions
    'load_model_and_tokenizer',
    'process_batches_on_gpu',
    'get_batch_embeddings',
    'calculate_single_sequence_embedding_last_hidden_state',
    'calculate_single_sequence_embedding_all_layers',
    'calculate_single_sequence_embedding_first_layer',
    'get_single_embedding_last_hidden_state',
    'get_single_embedding_all_layers',
    'get_single_embedding_first_layer',
] 