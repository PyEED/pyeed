import gc
import os
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
from transformers import EsmModel, EsmTokenizer

from pyeed.dbconnect import DatabaseConnector


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


def process_batches_on_gpu(
    data: list[tuple[str, str]],
    batch_size: int,
    model: Module,
    tokenizer: EsmTokenizer,
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
    logger.debug(f"Processing {len(data)} sequences on {device}.")

    model = model.to(device)

    # Split data into smaller batches
    for batch_start in range(0, len(data), batch_size):
        batch_end = min(batch_start + batch_size, len(data))
        batch = data[batch_start:batch_end]

        accessions, sequences = zip(*batch)

        current_batch_size = len(sequences)

        while current_batch_size > 0:
            try:
                # Compute embeddings
                embeddings_batch = get_batch_embeddings(
                    list(sequences[:current_batch_size]), model, tokenizer, device
                )

                # Update the database
                update_protein_embeddings_in_db(
                    db, list(accessions[:current_batch_size]), embeddings_batch
                )

                # Move to the next batch
                break  # Successful execution, move to the next batch

            except torch.cuda.OutOfMemoryError:
                torch.cuda.empty_cache()
                current_batch_size = max(
                    1, current_batch_size // 2
                )  # Reduce batch size
                logger.warning(
                    f"Reduced batch size to {current_batch_size} due to OOM error."
                )

    # Free memory
    del model
    torch.cuda.empty_cache()


def load_model_and_tokenizer(
    model_name: str,
    device: torch.device,
) -> Tuple[Any, Union[Any, None], torch.device]:
    """
    Loads the model and assigns it to a specific GPU.

    Args:
        model_name (str): The model name.
        device (str): The specific GPU device.

    Returns:
        Tuple: (model, tokenizer, device)
    """
    token = get_hf_token()
    tokenizer = None

    if "esmc" in model_name.lower():
        model = ESMC.from_pretrained(model_name)
    elif "esm3-sm-open-v1" in model_name.lower():
        model = ESM3.from_pretrained("esm3_sm_open_v1")
    else:
        full_model_name = (
            model_name
            if model_name.startswith("facebook/")
            else f"facebook/{model_name}"
        )
        model = EsmModel.from_pretrained(full_model_name, use_auth_token=token)
        tokenizer = EsmTokenizer.from_pretrained(full_model_name, use_auth_token=token)

    model = model.to(device)
    return model, tokenizer, device


def get_batch_embeddings(
    batch_sequences: list[str],
    model: Union[
        EsmModel,
        ESMC,
        DataParallel[Module],
        ESM3InferenceClient,
        ESM3,
    ],
    tokenizer_or_alphabet: Union[EsmTokenizer, None],
    device: torch.device,
    pool_embeddings: bool = True,
) -> list[NDArray[np.float64]]:
    """
    Generates mean-pooled embeddings for a batch of sequences.
    Supports ESM++, ESM-2 and ESM-3 models.

    Args:
        batch_sequences (list[str]): List of sequence strings.
        model: Loaded model (could be wrapped in DataParallel).
        tokenizer_or_alphabet: Tokenizer if needed.
        device: Inference device (CPU/GPU).
        pool_embeddings (bool): Whether to average embeddings across the sequence length.

    Returns:
        List of embeddings as NumPy arrays.
    """
    # First, determine the base model type
    base_model = model.module if isinstance(model, torch.nn.DataParallel) else model

    if isinstance(base_model, ESMC):
        # For ESMC models
        embedding_list = []
        with torch.no_grad():
            for sequence in batch_sequences:
                protein = ESMProtein(sequence=sequence)
                # Use the model directly - DataParallel handles internal distribution
                protein_tensor = base_model.encode(protein)
                logits_output = base_model.logits(
                    protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
                )
                if logits_output.embeddings is None:
                    raise ValueError(
                        "Model did not return embeddings. Check LogitsConfig settings."
                    )
                embeddings = logits_output.embeddings.cpu().numpy()
                if pool_embeddings:
                    embeddings = embeddings.mean(axis=1)
                embedding_list.append(embeddings[0])
        return embedding_list
    elif isinstance(base_model, ESM3):
        # For ESM3 models
        embedding_list = []
        with torch.no_grad():
            for sequence in batch_sequences:
                protein = ESMProtein(sequence=sequence)
                sequence_encoding = base_model.encode(protein)
                result = base_model.forward_and_sample(
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
    else:
        # ESM-2 logic
        assert tokenizer_or_alphabet is not None, "Tokenizer required for ESM-2 models"
        inputs = tokenizer_or_alphabet(
            batch_sequences, padding=True, truncation=True, return_tensors="pt"
        ).to(device)
        with torch.no_grad():
            outputs = model(**inputs, output_hidden_states=True)

        # Get last hidden state for each sequence
        hidden_states = outputs.last_hidden_state.cpu().numpy()

        if pool_embeddings:
            # Mean pooling across sequence length
            return [embedding.mean(axis=0) for embedding in hidden_states]
        return list(hidden_states)


def calculate_single_sequence_embedding_last_hidden_state(
    sequence: str,
    device: torch.device,
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
    model, tokenizer, device = load_model_and_tokenizer(model_name, device)
    return get_single_embedding_last_hidden_state(sequence, model, tokenizer, device)


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
    model, tokenizer, device = load_model_and_tokenizer(model_name, device)
    return get_single_embedding_all_layers(sequence, model, tokenizer, device)


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
    from esm.models.esmc import ESMC

    with torch.no_grad():
        if isinstance(model, ESMC):
            # ESM-3 logic
            from esm.sdk.api import ESMProtein, LogitsConfig

            protein = ESMProtein(sequence=sequence)
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

            embedding = (
                logits_output.hidden_states[-1][0].to(torch.float32).cpu().numpy()
            )
        else:
            # ESM-2 logic
            inputs = tokenizer(sequence, return_tensors="pt").to(device)
            outputs = model(**inputs)
            embedding = outputs.last_hidden_state[0, 1:-1, :].detach().cpu().numpy()

    # normalize the embedding
    embedding = embedding / np.linalg.norm(embedding, axis=1, keepdims=True)

    return embedding  # type: ignore


def get_single_embedding_all_layers(
    sequence: str, model: Any, tokenizer: Any, device: torch.device
) -> NDArray[np.float64]:
    """
    Generates normalized embeddings for each token in the sequence across all layers.

    For ESM-3 (ESMC) models, it assumes that passing
    LogitsConfig(return_hidden_states=True) returns a collection of layer embeddings.
    For ESM-2 models, it sets output_hidden_states=True.

    Args:
        sequence (str): The protein sequence to embed.
        model (Any): The transformer model to use.
        tokenizer (Any): The tokenizer for the model (None for ESMC).
        device (torch.device): The device to run the model on (CPU/GPU).

    Returns:
        NDArray[np.float64]: A numpy array containing the normalized token embeddings
        concatenated across all layers.
    """
    embeddings_list = []
    with torch.no_grad():
        if isinstance(model, ESMC):
            # For ESM-3: Use ESMProtein and request hidden states via LogitsConfig
            protein = ESMProtein(sequence=sequence)
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
                emb = emb / np.linalg.norm(emb, axis=1, keepdims=True)
                embeddings_list.append(emb)

        else:
            # For ESM-2: Get hidden states with output_hidden_states=True
            inputs = tokenizer(sequence, return_tensors="pt").to(device)
            outputs = model(**inputs, output_hidden_states=True)
            hidden_states = (
                outputs.hidden_states
            )  # Tuple: (layer0, layer1, ..., layerN)
            for layer_tensor in hidden_states:
                # Remove batch dimension and special tokens ([CLS] and [SEP])
                emb = layer_tensor[0, 1:-1, :].detach().cpu().numpy()
                emb = emb / np.linalg.norm(emb, axis=1, keepdims=True)
                embeddings_list.append(emb)

    return np.array(embeddings_list)


# The rest of your existing functions will need to be adapted in a similar way
# if they interact with the model or tokenizer directly


def free_memory() -> None:
    """
    Frees up memory by invoking garbage collection and clearing GPU caches.
    """
    gc.collect()
    if torch.backends.mps.is_available():
        torch.mps.empty_cache()
    elif torch.cuda.is_available():
        torch.cuda.empty_cache()


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
    # Prepare the data for batch update
    updates = [
        {"accession": acc, "embedding": emb.tolist()}
        for acc, emb in zip(accessions, embeddings_batch)
    ]

    # Cypher query for batch update
    query = """
    UNWIND $updates AS update
    MATCH (p:Protein {accession_id: update.accession})
    SET p.embedding = update.embedding
    """

    # Execute the update query with parameters
    db.execute_write(query, {"updates": updates})