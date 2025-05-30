"""
Utility functions for embedding operations.

Contains helper functions for token management, memory management, 
and sequence preprocessing.
"""

import gc
import os
import re

import torch
from huggingface_hub import HfFolder, login


def get_hf_token() -> str:
    """Get or request Hugging Face token."""
    if os.getenv("PYTEST_DISABLE_HF_LOGIN"):  # Disable Hugging Face login in tests
        return "dummy_token_for_tests"

    hf_folder = HfFolder()
    token = hf_folder.get_token()
    if not token:
        login()  # Login returns None, get token after login
        token = hf_folder.get_token()

    if isinstance(token, str):
        return token
    else:
        raise RuntimeError("Failed to get Hugging Face token")


def preprocess_sequence_for_prott5(sequence: str) -> str:
    """
    Preprocesses a protein sequence for ProtT5 models.
    
    Args:
        sequence: Raw protein sequence
        
    Returns:
        Preprocessed sequence with spaces between amino acids and rare AAs mapped to X
    """
    # Map rare amino acids to X and add spaces between amino acids
    sequence = re.sub(r"[UZOB]", "X", sequence.upper())
    return " ".join(list(sequence))


def free_memory() -> None:
    """
    Frees up memory by invoking garbage collection and clearing GPU caches.
    """
    gc.collect()
    if torch.backends.mps.is_available():
        torch.mps.empty_cache()
    elif torch.cuda.is_available():
        torch.cuda.empty_cache()


def determine_model_type(model_name: str) -> str:
    """
    Determine the model type based on model name.
    
    Args:
        model_name: Name of the model
        
    Returns:
        Model type string
    """
    model_name_lower = model_name.lower()
    
    if "esmc" in model_name_lower:
        return "esmc"
    elif "esm3" in model_name_lower:
        return "esm3"
    elif "prot_t5" in model_name_lower or "prott5" in model_name_lower:
        return "prott5"
    else:
        return "esm2"  # Default to ESM-2 for other facebook/esm models 