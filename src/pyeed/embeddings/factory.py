"""
Factory for creating embedding model instances.

Provides a centralized way to create different types of embedding models
based on model names and automatically handles device assignment.
"""

from typing import Any, Tuple, Union

import torch
from torch.nn import DataParallel, Module

from .base import BaseEmbeddingModel
from .models import (
    ESM2EmbeddingModel,
    ESM3EmbeddingModel,
    ESMCEmbeddingModel,
    ProtT5EmbeddingModel,
)
from .utils import determine_model_type


class ModelFactory:
    """Factory for creating embedding model instances."""
    
    @staticmethod
    def create_model(
        model_name: str, 
        device: torch.device = torch.device("cuda:0")
    ) -> BaseEmbeddingModel:
        """
        Create an embedding model instance based on the model name.
        
        Args:
            model_name: Name of the model to create
            device: Device to run the model on
            
        Returns:
            BaseEmbeddingModel instance
        """
        model_type = determine_model_type(model_name)
        
        if model_type == "esmc":
            return ESMCEmbeddingModel(model_name, device)
        elif model_type == "esm3":
            return ESM3EmbeddingModel(model_name, device)
        elif model_type == "prott5":
            return ProtT5EmbeddingModel(model_name, device)
        else:  # Default to ESM-2
            return ESM2EmbeddingModel(model_name, device)
    
    @staticmethod
    def load_model_and_tokenizer(
        model_name: str,
        device: torch.device = torch.device("cuda:0"),
    ) -> Tuple[Union[Any, DataParallel[Module]], Union[Any, None], torch.device]:
        """
        Load model and tokenizer using the factory pattern.
        
        This method maintains compatibility with the original function signature
        while using the new OOP structure internally.
        
        Args:
            model_name: The model name
            device: The specific GPU device
            
        Returns:
            Tuple: (model, tokenizer, device)
        """
        embedding_model = ModelFactory.create_model(model_name, device)
        model, tokenizer = embedding_model.load_model()
        
        return model, tokenizer, device 