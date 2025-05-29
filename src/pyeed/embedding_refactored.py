"""
Refactored embedding module that maintains original function signatures.

This module provides the same interface as the original embedding.py while
using the new organized structure with model classes, factory, and processor.
"""

import gc
import os
import re
from typing import Any, Tuple, Union

import numpy as np
import torch
from esm.models.esm3 import ESM3
from esm.models.esmc import ESMC
from esm.sdk.api import ESM3InferenceClient, ESMProtein, LogitsConfig, SamplingConfig
from huggingface_hub import HfFolder, login
from loguru import logger
from numpy.typing import NDArray
from torch.nn import DataParallel, Module
from transformers import EsmModel, EsmTokenizer, T5Model, T5Tokenizer

from pyeed.dbconnect import DatabaseConnector
from pyeed.embeddings.processor import get_processor
from pyeed.embeddings.factory import ModelFactory
from pyeed.embeddings.database import update_protein_embeddings_in_db as _update_protein_embeddings_in_db
from pyeed.embeddings.utils import get_hf_token as _get_hf_token, preprocess_sequence_for_prott5 as _preprocess_sequence_for_prott5, free_memory as _free_memory


# ============================================================================
# Original function signatures maintained for backward compatibility
# ============================================================================

def get_hf_token() -> str:
    """Get or request Hugging Face token."""
    return _get_hf_token()


def process_batches_on_gpu(
    data: list[tuple[str, str]],
    batch_size: int,
    model: Union[EsmModel, ESMC, ESM3, T5Model, DataParallel[Module]],
    tokenizer: Union[EsmTokenizer, T5Tokenizer, None],
    db: DatabaseConnector,
    device: torch.device,
) -> None:
    """
    Splits data into batches and processes them on a single GPU.

    Args:
        data (list): List of (accession_id, sequence) tuples.
        batch_size (int): Size of each batch.
        model: The model instance for this GPU.
        tokenizer: The tokenizer for the model.
        device (str): The assigned GPU device.
        db: Database connection.
    """
    processor = get_processor()
    processor.process_batches_on_gpu(data, batch_size, model, tokenizer, db, device)


def load_model_and_tokenizer(
    model_name: str,
    device: torch.device = torch.device("cuda:0"),
) -> Tuple[Union[EsmModel, ESMC, ESM3, T5Model], Union[EsmTokenizer, T5Tokenizer, None], torch.device]:
    """
    Loads the model and assigns it to a specific GPU.

    Args:
        model_name (str): The model name.
        device (str): The specific GPU device.

    Returns:
        Tuple: (model, tokenizer, device)
    """
    return ModelFactory.load_model_and_tokenizer(model_name, device)


def preprocess_sequence_for_prott5(sequence: str) -> str:
    """
    Preprocesses a protein sequence for ProtT5 models.
    
    Args:
        sequence: Raw protein sequence
        
    Returns:
        Preprocessed sequence with spaces between amino acids and rare AAs mapped to X
    """
    return _preprocess_sequence_for_prott5(sequence)


def get_batch_embeddings(
    batch_sequences: list[str],
    model: Union[
        EsmModel,
        ESMC,
        DataParallel[Module],
        ESM3,
        T5Model,
    ],
    tokenizer_or_alphabet: Union[EsmTokenizer, T5Tokenizer, None],
    device: torch.device,
    pool_embeddings: bool = True,
) -> list[NDArray[np.float64]]:
    """
    Generates mean-pooled embeddings for a batch of sequences.
    Supports ESM++, ESM-2, ESM-3 and ProtT5 models.

    Args:
        batch_sequences (list[str]): List of sequence strings.
        model: Loaded model (could be wrapped in DataParallel).
        tokenizer_or_alphabet: Tokenizer if needed.
        device: Inference device (CPU/GPU).
        pool_embeddings (bool): Whether to average embeddings across the sequence length.

    Returns:
        List of embeddings as NumPy arrays.
    """
    processor = get_processor()
    return processor.get_batch_embeddings_unified(
        batch_sequences, model, tokenizer_or_alphabet, device, pool_embeddings
    )


def calculate_single_sequence_embedding_last_hidden_state(
    sequence: str,
    device: torch.device = torch.device("cuda:0"),
    model_name: str = "facebook/esm2_t33_650M_UR50D",
) -> NDArray[np.float64]:
    """
    Calculates an embedding for a single sequence.

    Args:
        sequence: Input protein sequence
        model_name: Name of the ESM model to use

    Returns:
        NDArray[np.float64]: Normalized embedding vector for the sequence
    """
    processor = get_processor()
    return processor.calculate_single_sequence_embedding_last_hidden_state(
        sequence, device, model_name
    )


def calculate_single_sequence_embedding_all_layers(
    sequence: str,
    device: torch.device,
    model_name: str = "facebook/esm2_t33_650M_UR50D",
) -> NDArray[np.float64]:
    """
    Calculates embeddings for a single sequence across all layers.

    Args:
        sequence: Input protein sequence
        model_name: Name of the ESM model to use

    Returns:
        NDArray[np.float64]: A numpy array containing layer embeddings for the sequence.
    """
    processor = get_processor()
    return processor.calculate_single_sequence_embedding_all_layers(
        sequence, device, model_name
    )


def get_single_embedding_last_hidden_state(
    sequence: str, model: Any, tokenizer: Any, device: torch.device
) -> NDArray[np.float64]:
    """Generate embeddings for a single sequence using the last hidden state.

    Args:
        sequence (str): The protein sequence to embed
        model (Any): The transformer model to use
        tokenizer (Any): The tokenizer for the model
        device (torch.device): The device to run the model on (CPU/GPU)

    Returns:
        np.ndarray: Normalized embeddings for each token in the sequence
    """
    processor = get_processor()
    return processor.get_single_embedding_last_hidden_state(sequence, model, tokenizer, device)


def get_single_embedding_all_layers(
    sequence: str, model: Any, tokenizer: Any, device: torch.device
) -> NDArray[np.float64]:
    """
    Generates normalized embeddings for each token in the sequence across all layers.

    For ESM-3 (ESMC) models, it assumes that passing
    LogitsConfig(return_hidden_states=True) returns a collection of layer embeddings.
    For ESM-2 models, it sets output_hidden_states=True.
    For ProtT5 models, it gets encoder hidden states.

    Args:
        sequence (str): The protein sequence to embed.
        model (Any): The transformer model to use.
        tokenizer (Any): The tokenizer for the model (None for ESMC).
        device (torch.device): The device to run the model on (CPU/GPU).

    Returns:
        NDArray[np.float64]: A numpy array containing the normalized token embeddings
        concatenated across all layers.
    """
    processor = get_processor()
    return processor.get_single_embedding_all_layers(sequence, model, tokenizer, device)


def calculate_single_sequence_embedding_first_layer(
    sequence: str, model_name: str = "facebook/esm2_t33_650M_UR50D", device: torch.device = torch.device("cuda:0"),
) -> NDArray[np.float64]:
    """
    Calculates an embedding for a single sequence using the first layer.
    """
    processor = get_processor()
    return processor.calculate_single_sequence_embedding_first_layer(sequence, model_name, device)


def get_single_embedding_first_layer(
    sequence: str, model: Any, tokenizer: Any, device: torch.device
) -> NDArray[np.float64]:
    """
    Generates normalized embeddings for each token in the sequence using the first layer.
    """
    processor = get_processor()
    return processor.get_single_embedding_first_layer(sequence, model, tokenizer, device)


def free_memory() -> None:
    """
    Frees up memory by invoking garbage collection and clearing GPU caches.
    """
    _free_memory()


def update_protein_embeddings_in_db(
    db: DatabaseConnector,
    accessions: list[str],
    embeddings_batch: list[NDArray[np.float64]],
) -> None:
    """
    Updates the embeddings for a batch of proteins in the database.

    Args:
        db (DatabaseConnector): The database connector.
        accessions (list[str]): The accessions of the proteins to update.
        embeddings_batch (list[NDArray[np.float64]]): The embeddings to update.
    """
    _update_protein_embeddings_in_db(db, accessions, embeddings_batch) 