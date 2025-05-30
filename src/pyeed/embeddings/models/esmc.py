"""
ESMC model implementation for protein embeddings.
"""

from typing import List, Tuple, Optional, cast
import torch
import numpy as np
from numpy.typing import NDArray
from loguru import logger
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig

from ..base import BaseEmbeddingModel, normalize_embedding


class ESMCEmbeddingModel(BaseEmbeddingModel):
    """ESMC model implementation."""
    
    def __init__(self, model_name: str, device: torch.device):
        super().__init__(model_name, device)
    
    def load_model(self) -> Tuple[ESMC, None]:
        """Load ESMC model with improved error handling."""
        try:
            # Try to disable tqdm to avoid threading issues
            import os
            os.environ['DISABLE_TQDM'] = 'True'
            
            model = ESMC.from_pretrained(self.model_name)
            model = model.to(self.device)
            
            self.model = model
            
            return model, None
            
        except Exception as e:
            if "tqdm" in str(e).lower() or "_lock" in str(e).lower():
                logger.warning(f"ESMC model loading failed due to tqdm threading issue: {e}. Retrying with threading workaround...")
                
                # Try alternative approach with threading lock
                import threading
                import time
                
                # Add a small delay and retry
                time.sleep(0.1 + torch.cuda.current_device() * 0.05)  # Staggered delay per GPU
                
                try:
                    # Try importing tqdm and resetting its state
                    try:
                        import tqdm
                        if hasattr(tqdm.tqdm, '_lock'):
                            delattr(tqdm.tqdm, '_lock')
                    except:
                        pass
                    
                    model = ESMC.from_pretrained(self.model_name)
                    model = model.to(self.device)
                    
                    self.model = model
                    
                    return model, None
                    
                except Exception as retry_error:
                    logger.error(f"ESMC model loading failed even after retry: {retry_error}")
                    raise retry_error
            else:
                logger.error(f"ESMC model loading failed: {e}")
                raise e
    
    def preprocess_sequence(self, sequence: str) -> ESMProtein:
        """ESMC uses ESMProtein objects."""
        return ESMProtein(sequence=sequence)
    
    def get_batch_embeddings(
        self, 
        sequences: List[str], 
        pool_embeddings: bool = True
    ) -> List[NDArray[np.float64]]:
        """Get embeddings for a batch of sequences using ESMC."""
        if self.model is None:
            self.load_model()
        
        # Type cast to ensure type checker knows it's not None
        model = cast(ESMC, self.model)
        
        embedding_list = []
        with torch.no_grad():
            for sequence in sequences:
                protein = self.preprocess_sequence(sequence)
                # Use the model directly - DataParallel handles internal distribution
                protein_tensor = model.encode(protein)
                logits_output = model.logits(
                    protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
                )
                if logits_output.embeddings is None:
                    raise ValueError(
                        "Model did not return embeddings. Check LogitsConfig settings."
                    )
                embeddings = logits_output.embeddings.cpu().numpy()
                # drop the special tokens
                embeddings = embeddings[:, 1:-1, :]
                if pool_embeddings:
                    embeddings = embeddings.mean(axis=1)
                embedding_list.append(embeddings[0])
        return embedding_list
    
    def get_single_embedding_last_hidden_state(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get last hidden state embedding for a single sequence."""
        if self.model is None:
            self.load_model()
        
        # Type cast to ensure type checker knows it's not None
        model = cast(ESMC, self.model)
        
        with torch.no_grad():
            protein = self.preprocess_sequence(sequence)
            protein_tensor = model.encode(protein)
            logits_output = model.logits(
                protein_tensor,
                LogitsConfig(
                    sequence=True,
                    return_embeddings=True,
                    return_hidden_states=True,
                ),
            )
            # Ensure hidden_states is not None before accessing it
            if logits_output.hidden_states is None:
                raise ValueError(
                    "Model did not return hidden states. Check LogitsConfig settings."
                )

            # remove special tokens
            embedding = (
                logits_output.hidden_states[-1][0][1:-1].to(torch.float32).cpu().numpy()
            )

        # Normalize the embedding
        embedding = normalize_embedding(embedding)
        return embedding
    
    def get_single_embedding_all_layers(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get embeddings from all layers for a single sequence."""
        if self.model is None:
            self.load_model()
        
        # Type cast to ensure type checker knows it's not None
        model = cast(ESMC, self.model)
        
        embeddings_list = []
        with torch.no_grad():
            protein = self.preprocess_sequence(sequence)
            protein_tensor = model.encode(protein)
            logits_output = model.logits(
                protein_tensor,
                LogitsConfig(
                    sequence=True,
                    return_embeddings=True,
                    return_hidden_states=True,
                ),
            )
            # Ensure hidden_states is not None before iterating
            if logits_output.hidden_states is None:
                raise ValueError(
                    "Model did not return hidden states. Check if return_hidden_states=True is supported."
                )

            # logits_output.hidden_states should be a tuple of tensors: (layer, batch, seq_len, hidden_dim)
            for layer_tensor in logits_output.hidden_states:
                # Remove batch dimension and (if applicable) any special tokens
                emb = layer_tensor[0].to(torch.float32).cpu().numpy()
                # If your model adds special tokens, adjust the slicing (e.g., emb[1:-1])
                emb = normalize_embedding(emb)
                embeddings_list.append(emb)

        return np.array(embeddings_list)
    
    def get_single_embedding_first_layer(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """Get first layer embedding for a single sequence."""
        if self.model is None:
            self.load_model()
        
        # Type cast to ensure type checker knows it's not None
        model = cast(ESMC, self.model)
        
        with torch.no_grad():
            protein = self.preprocess_sequence(sequence)
            protein_tensor = model.encode(protein)
            logits_output = model.logits(
                protein_tensor,
                LogitsConfig(
                    sequence=True,
                    return_embeddings=True,
                    return_hidden_states=True,
                ),
            )
            if logits_output.hidden_states is None:
                raise ValueError(
                    "Model did not return hidden states. Check LogitsConfig settings."
                )
            embedding = (
                logits_output.hidden_states[0][0].to(torch.float32).cpu().numpy()
            )

        # Normalize the embedding
        embedding = normalize_embedding(embedding)
        return embedding 
    
    def get_final_embeddings(
        self, 
        sequence: str
    ) -> NDArray[np.float64]:
        """
        Get final embeddings for ESMC with robust fallback.
        
        Provides a more robust embedding extraction that prioritizes
        batch embeddings (properly pooled) over last hidden state.
        """
        try:
            # For ESMC, batch embeddings with pooling is more reliable and memory efficient
            embeddings = self.get_batch_embeddings([sequence], pool_embeddings=True)
            if embeddings and len(embeddings) > 0:
                return np.asarray(embeddings[0], dtype=np.float64)
            else:
                raise ValueError("Batch embeddings method returned empty results")
        except (torch.cuda.OutOfMemoryError, RuntimeError) as e:
            if "out of memory" in str(e).lower():
                logger.warning(f"Batch embeddings method failed due to OOM for ESMC: {e}. Clearing cache and trying minimal approach.")
                # Clear cache and try a more memory-efficient approach
                torch.cuda.empty_cache()
                try:
                    # Minimal approach - just get embeddings without requesting hidden states
                    if self.model is None:
                        self.load_model()
                    
                    model = cast(ESMC, self.model)
                    
                    with torch.no_grad():
                        protein = self.preprocess_sequence(sequence)
                        protein_tensor = model.encode(protein)
                        logits_output = model.logits(
                            protein_tensor, 
                            LogitsConfig(sequence=True, return_embeddings=True)
                        )
                        if logits_output.embeddings is None:
                            raise ValueError("Model did not return embeddings")
                        
                        # Get embeddings and pool them properly
                        embeddings = logits_output.embeddings.cpu().numpy()
                        logger.info(f"Embeddings shape: {embeddings.shape}")
                        
                        # Pool across sequence dimension to get single vector
                        pooled_embedding = embeddings.mean(axis=1)[0]
                        
                        return np.asarray(pooled_embedding, dtype=np.float64)
                        
                except Exception as minimal_error:
                    logger.error(f"Minimal embedding extraction also failed for ESMC: {minimal_error}")
                    raise ValueError(f"ESMC embedding extraction failed with OOM: {minimal_error}")
            else:
                raise e
        except Exception as e:
            logger.error(f"All embedding extraction methods failed for ESMC: {e}")
            raise ValueError(f"ESMC embedding extraction failed: {e}") 