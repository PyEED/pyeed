"""
Organized embedding module for protein language models.

This module provides both the new organized structure and backward compatibility
with the original embedding.py interface.
"""

# New organized structure
from .base import BaseEmbeddingModel, ModelType, normalize_embedding
from .factory import ModelFactory
from .processor import EmbeddingProcessor, get_processor
from .utils import get_hf_token, preprocess_sequence_for_prott5, free_memory, determine_model_type
from .database import update_protein_embeddings_in_db
from .models import ESM2EmbeddingModel, ESMCEmbeddingModel, ESM3EmbeddingModel, ProtT5EmbeddingModel

# Backward compatibility imports from old embedding.py
try:
    from ..embedding import (
        load_model_and_tokenizer,
        process_batches_on_gpu,
        get_batch_embeddings,
        calculate_single_sequence_embedding_last_hidden_state,
        calculate_single_sequence_embedding_all_layers,
        calculate_single_sequence_embedding_first_layer,
        get_single_embedding_last_hidden_state,
        get_single_embedding_all_layers,
        get_single_embedding_first_layer
    )
except ImportError:
    # If old embedding.py is not available, use processor methods for compatibility
    _processor = get_processor()
    
    def load_model_and_tokenizer(model_name: str, device=None):
        """Backward compatibility function."""
        # This is handled internally by the processor now
        return None, None, device
    
    def process_batches_on_gpu(data, batch_size, model, tokenizer, db, device):
        """Backward compatibility function."""
        return _processor.process_batches_on_gpu(data, batch_size, model, tokenizer, db, device)
    
    def get_batch_embeddings(batch_sequences, model, tokenizer, device, pool_embeddings=True):
        """Backward compatibility function."""
        return _processor.get_batch_embeddings_unified(batch_sequences, model, tokenizer, device, pool_embeddings)
    
    def calculate_single_sequence_embedding_last_hidden_state(sequence, device=None, model_name="facebook/esm2_t33_650M_UR50D"):
        """Backward compatibility function."""
        return _processor.calculate_single_embedding(sequence, model_name, "last_hidden_state", device)
    
    def calculate_single_sequence_embedding_all_layers(sequence, device, model_name="facebook/esm2_t33_650M_UR50D"):
        """Backward compatibility function."""
        return _processor.calculate_single_embedding(sequence, model_name, "all_layers", device)
    
    def calculate_single_sequence_embedding_first_layer(sequence, model_name="facebook/esm2_t33_650M_UR50D", device=None):
        """Backward compatibility function."""
        return _processor.calculate_single_embedding(sequence, model_name, "first_layer", device)
    
    def get_single_embedding_last_hidden_state(sequence, model, tokenizer, device):
        """Backward compatibility function."""
        return _processor.get_single_embedding_last_hidden_state(sequence, model, tokenizer, device)
    
    def get_single_embedding_all_layers(sequence, model, tokenizer, device):
        """Backward compatibility function."""
        return _processor.get_single_embedding_all_layers(sequence, model, tokenizer, device)
    
    def get_single_embedding_first_layer(sequence, model, tokenizer, device):
        """Backward compatibility function."""
        return _processor.get_single_embedding_first_layer(sequence, model, tokenizer, device)

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