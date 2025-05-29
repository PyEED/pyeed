"""
ESM-2 model implementation for protein embeddings.
"""

from typing import List, Tuple, Optional, Any, cast
import torch
import numpy as np
from numpy.typing import NDArray
from transformers import EsmModel, EsmTokenizer
from loguru import logger

from ..base import BaseEmbeddingModel, normalize_embedding
from ..utils import get_hf_token


class ESM2EmbeddingModel(BaseEmbeddingModel):
    """ESM-2 model implementation."""
    
    def __init__(self, model_name: str, device: torch.device):
        super().__init__(model_name, device)
    
    def load_model(self) -> Tuple[EsmModel, EsmTokenizer]:
        """Load ESM-2 model and tokenizer."""
        token = get_hf_token()
        
        full_model_name = (
            self.model_name
            if self.model_name.startswith("facebook/")
            else f"facebook/{self.model_name}"
        )
        
        model = EsmModel.from_pretrained(full_model_name, use_auth_token=token)
        tokenizer = EsmTokenizer.from_pretrained(full_model_name, use_auth_token=token)
        
        # Move to device
        model = model.to(self.device)
        
        self.model = model
        self.tokenizer = tokenizer
        
        return model, tokenizer
    
    def preprocess_sequence(self, sequence: str) -> str:
        """ESM-2 doesn't need special preprocessing."""
        return sequence
    
    def get_batch_embeddings(
        self, 
        sequences: List[str], 
        pool_embeddings: bool = True
    ) -> List[NDArray[np.float64]]:
        """Get embeddings for a batch of sequences using ESM-2."""
        if self.model is None or self.tokenizer is None:
            self.load_model()
        
        # Type cast to ensure type checker knows they're not None
        model = cast(EsmModel, self.model)
        tokenizer = cast(EsmTokenizer, self.tokenizer)
        
        inputs = tokenizer(
            sequences, padding=True, truncation=True, return_tensors="pt"
        ).to(self.device)
        
        with torch.no_grad():
            outputs = model(**inputs, output_hidden_states=True)

        # Get last hidden state for each sequence
        hidden_states = outputs.last_hidden_state.cpu().numpy()

        if pool_embeddings:
            # Mean pooling across sequence length
            return [embedding.mean(axis=0) for embedding in hidden_states]
        return list(hidden_states)
    
    def get_single_embedding_last_hidden_state(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get last hidden state embedding for a single sequence."""
        if self.model is None or self.tokenizer is None:
            self.load_model()
        
        # Type cast to ensure type checker knows they're not None
        model = cast(EsmModel, self.model)
        tokenizer = cast(EsmTokenizer, self.tokenizer)
        
        inputs = tokenizer(sequence, return_tensors="pt").to(self.device)
        
        with torch.no_grad():
            outputs = model(**inputs)
        
        # Remove batch dimension and special tokens ([CLS] and [SEP])
        embedding = outputs.last_hidden_state[0, 1:-1, :].detach().cpu().numpy()
        return embedding
    
    def get_single_embedding_all_layers(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get embeddings from all layers for a single sequence."""
        if self.model is None or self.tokenizer is None:
            self.load_model()
        
        # Type cast to ensure type checker knows they're not None
        model = cast(EsmModel, self.model)
        tokenizer = cast(EsmTokenizer, self.tokenizer)
        
        inputs = tokenizer(sequence, return_tensors="pt").to(self.device)
        
        with torch.no_grad():
            outputs = model(**inputs, output_hidden_states=True)
        
        embeddings_list = []
        hidden_states = outputs.hidden_states  # Tuple: (layer0, layer1, ..., layerN)
        
        for layer_tensor in hidden_states:
            # Remove batch dimension and special tokens ([CLS] and [SEP])
            emb = layer_tensor[0, 1:-1, :].detach().cpu().numpy()
            emb = normalize_embedding(emb)
            embeddings_list.append(emb)

        return np.array(embeddings_list)
    
    def get_single_embedding_first_layer(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get first layer embedding for a single sequence."""
        if self.model is None or self.tokenizer is None:
            self.load_model()
        
        # Type cast to ensure type checker knows they're not None
        model = cast(EsmModel, self.model)
        tokenizer = cast(EsmTokenizer, self.tokenizer)
        
        inputs = tokenizer(sequence, return_tensors="pt").to(self.device)
        
        with torch.no_grad():
            outputs = model(**inputs, output_hidden_states=True)
        
        # Get the first layer's hidden states for all residues (excluding special tokens)
        embedding = outputs.hidden_states[0][0, 1:-1, :].detach().cpu().numpy()
        
        # Normalize the embedding
        embedding = normalize_embedding(embedding)
        return embedding

    def get_final_embeddings(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """
        Get final embeddings for ESM-2 with robust fallback.
        
        Provides a more robust embedding extraction that prioritizes
        batch processing for better performance.
        """
        try:
            # For ESM-2, batch processing is more efficient
            embeddings = self.get_batch_embeddings([sequence], pool_embeddings=True)
            if embeddings and len(embeddings) > 0:
                return embeddings[0]
            else:
                raise ValueError("Batch embeddings method returned empty results")
        except Exception as e:
            logger.warning(f"Batch embeddings method failed for ESM-2: {e}. Trying single sequence method.")
            try:
                # Fallback to single sequence method
                return self.get_single_embedding_last_hidden_state(sequence)
            except Exception as fallback_error:
                logger.error(f"All embedding extraction methods failed for ESM-2: {fallback_error}")
                raise ValueError(f"ESM-2 embedding extraction failed: {fallback_error}") 