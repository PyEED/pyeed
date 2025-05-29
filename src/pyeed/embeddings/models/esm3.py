"""
ESM-3 model implementation for protein embeddings.
"""

from typing import List, Tuple, Optional, cast
import torch
import numpy as np
from numpy.typing import NDArray
from loguru import logger
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein, SamplingConfig

from ..base import BaseEmbeddingModel, normalize_embedding


class ESM3EmbeddingModel(BaseEmbeddingModel):
    """ESM-3 model implementation."""
    
    def __init__(self, model_name: str, device: torch.device):
        super().__init__(model_name, device)
    
    def load_model(self) -> Tuple[ESM3, None]:
        """Load ESM3 model."""
        model = ESM3.from_pretrained("esm3_sm_open_v1")
        model = model.to(self.device)
        
        self.model = model
        
        return model, None
    
    def preprocess_sequence(self, sequence: str) -> ESMProtein:
        """ESM3 uses ESMProtein objects."""
        return ESMProtein(sequence=sequence)
    
    def get_batch_embeddings(
        self, 
        sequences: List[str], 
        pool_embeddings: bool = True
    ) -> List[NDArray[np.float64]]:
        """Get embeddings for a batch of sequences using ESM3."""
        if self.model is None:
            self.load_model()
        
        # Type cast to ensure type checker knows it's not None
        model = cast(ESM3, self.model)
        
        embedding_list = []
        with torch.no_grad():
            for sequence in sequences:
                protein = self.preprocess_sequence(sequence)
                sequence_encoding = model.encode(protein)
                result = model.forward_and_sample(
                    sequence_encoding,
                    SamplingConfig(return_per_residue_embeddings=True),
                )
                if result is None or result.per_residue_embedding is None:
                    raise ValueError("Model did not return embeddings")
                embeddings = (
                    result.per_residue_embedding.to(torch.float32).cpu().numpy()
                )
                if pool_embeddings:
                    embeddings = embeddings.mean(axis=0)
                embedding_list.append(embeddings)
        return embedding_list
    
    def get_single_embedding_last_hidden_state(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get last hidden state embedding for a single sequence."""
        if self.model is None:
            self.load_model()
        
        # Type cast to ensure type checker knows it's not None
        model = cast(ESM3, self.model)
        
        with torch.no_grad():
            protein = self.preprocess_sequence(sequence)
            sequence_encoding = model.encode(protein)
            embedding = model.forward_and_sample(
                sequence_encoding,
                SamplingConfig(return_per_residue_embeddings=True),
            )
            if embedding is None or embedding.per_residue_embedding is None:
                raise ValueError("Model did not return embeddings")
            embedding = embedding.per_residue_embedding.to(torch.float32).cpu().numpy()

        # Normalize the embedding
        embedding = normalize_embedding(embedding)
        return embedding
    
    def get_single_embedding_all_layers(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get embeddings from all layers for a single sequence."""
        # ESM3 doesn't support all layers extraction in the same way
        # This is a simplified implementation - might need enhancement based on ESM3 capabilities
        if self.model is None:
            self.load_model()
        
        # Type cast to ensure type checker knows it's not None
        model = cast(ESM3, self.model)
        
        with torch.no_grad():
            protein = self.preprocess_sequence(sequence)
            sequence_encoding = model.encode(protein)
            result = model.forward_and_sample(
                sequence_encoding,
                SamplingConfig(return_per_residue_embeddings=True),
            )
            if result is None or result.per_residue_embedding is None:
                raise ValueError("Model did not return embeddings")
            
            # For ESM3, we return the per-residue embedding as a single layer
            # This might need adjustment based on actual ESM3 API capabilities
            embedding = result.per_residue_embedding.to(torch.float32).cpu().numpy()
            embedding = normalize_embedding(embedding)

        # Return as a single layer array for consistency with other models
        return np.array([embedding])
    
    def get_single_embedding_first_layer(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get first layer embedding for a single sequence."""
        # For ESM3, this is the same as the per-residue embedding
        if self.model is None:
            self.load_model()
        
        # Type cast to ensure type checker knows it's not None
        model = cast(ESM3, self.model)
        
        with torch.no_grad():
            protein = self.preprocess_sequence(sequence)
            sequence_encoding = model.encode(protein)
            result = model.forward_and_sample(
                sequence_encoding,
                SamplingConfig(return_per_residue_embeddings=True),
            )
            if result is None or result.per_residue_embedding is None:
                raise ValueError("Model did not return embeddings")
            embedding = result.per_residue_embedding.to(torch.float32).cpu().numpy()

        # Normalize the embedding
        embedding = normalize_embedding(embedding)
        return embedding 
    
    def get_final_embeddings(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """
        Get final embeddings for ESM3 with robust fallback.
        
        ESM3 has different API structure, so this provides a more robust
        embedding extraction that works reliably across different ESM3 versions.
        """
        try:
            # Try to get the standard per-residue embedding
            return self.get_single_embedding_last_hidden_state(sequence)
        except Exception as e:
            # If that fails, try alternative method
            logger.warning(f"Standard embedding method failed for ESM3: {e}. Trying alternative method.")
            try:
                if self.model is None:
                    self.load_model()
                
                model = cast(ESM3, self.model)
                
                with torch.no_grad():
                    protein = self.preprocess_sequence(sequence)
                    sequence_encoding = model.encode(protein)
                    # Try with minimal sampling config
                    result = model.forward_and_sample(
                        sequence_encoding,
                        SamplingConfig()
                    )
                    
                    # Extract any available embedding
                    if hasattr(result, 'per_residue_embedding') and result.per_residue_embedding is not None:
                        embedding = result.per_residue_embedding.to(torch.float32).cpu().numpy()
                        return embedding
                    else:
                        # Last resort: use a simple mean-pooled sequence representation
                        logger.warning("No per-residue embeddings available, using basic fallback")
                        raise ValueError("Could not extract any embeddings from ESM3 model")
            except Exception as fallback_error:
                logger.error(f"All embedding extraction methods failed for ESM3: {fallback_error}")
                raise ValueError(f"ESM3 embedding extraction failed: {fallback_error}") 