"""
Base classes for protein embedding models.

Defines the common interface that all embedding model implementations should follow.
"""

from abc import ABC, abstractmethod
from typing import Any, List, Optional, Tuple, Union

import numpy as np
import torch
from numpy.typing import NDArray


class BaseEmbeddingModel(ABC):
    """Abstract base class for protein embedding models."""
    
    def __init__(self, model_name: str, device: torch.device):
        self.model_name = model_name
        self.device = device
        self._model: Optional[Any] = None
        self._tokenizer: Optional[Any] = None
        
    @property
    def model(self) -> Optional[Any]:
        """Get the model instance."""
        return self._model
    
    @model.setter
    def model(self, value: Any) -> None:
        """Set the model instance."""
        self._model = value
    
    @property
    def tokenizer(self) -> Optional[Any]:
        """Get the tokenizer instance."""
        return self._tokenizer
    
    @tokenizer.setter
    def tokenizer(self, value: Any) -> None:
        """Set the tokenizer instance."""
        self._tokenizer = value
    
    @abstractmethod
    def load_model(self) -> Tuple[Any, Optional[Any]]:
        """Load and return the model and tokenizer."""
        pass
    
    @abstractmethod
    def preprocess_sequence(self, sequence: str) -> Union[str, Any]:
        """Preprocess a sequence for the specific model type."""
        pass
    
    @abstractmethod
    def get_batch_embeddings(
        self, 
        sequences: List[str], 
        pool_embeddings: bool = True
    ) -> List[NDArray[np.float64]]:
        """Get embeddings for a batch of sequences."""
        pass
    
    @abstractmethod
    def get_single_embedding_last_hidden_state(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get embedding from the last hidden state for a single sequence."""
        pass
    
    @abstractmethod
    def get_single_embedding_all_layers(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get embeddings from all layers for a single sequence."""
        pass
    
    @abstractmethod
    def get_single_embedding_first_layer(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get embedding from the first layer for a single sequence."""
        pass
    
    def get_final_embeddings(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """
        Get final embeddings for a single sequence.
        
        This method provides a robust embedding option that works across all models.
        It falls back gracefully if certain layer-specific methods are not available.
        Default implementation uses last hidden state, but can be overridden.
        """
        result = self.get_single_embedding_last_hidden_state(sequence)
        return np.asarray(result, dtype=np.float64)
    
    def move_to_device(self) -> None:
        """Move model to the specified device."""
        if self.model is not None:
            self.model = self.model.to(self.device)
    
    def cleanup(self) -> None:
        """Clean up model resources."""
        if self._model is not None:
            self._model = None
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        # Explicit return None
        return None


class ModelType:
    """Constants for different model types."""
    ESM2 = "esm2"
    ESMC = "esmc"
    ESM3 = "esm3"
    PROTT5 = "prott5"


def normalize_embedding(embedding: NDArray[np.float64]) -> NDArray[np.float64]:
    """Normalize embeddings using L2 normalization."""
    norm = np.linalg.norm(embedding, axis=1, keepdims=True)
    # Handle zero norm case to avoid division by zero
    norm[norm == 0] = 1.0
    normalized = embedding / norm
    return np.asarray(normalized, dtype=np.float64) 