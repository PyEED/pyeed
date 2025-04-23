import gc
import logging
import os
from typing import Any, Tuple, Union

import numpy as np
import torch
from esm.models.esm3 import ESM3
from esm.models.esmc import ESMC
from esm.sdk.api import ESM3InferenceClient, ESMProtein, LogitsConfig, SamplingConfig
from huggingface_hub import HfFolder, login
from numpy.typing import NDArray
from transformers import EsmModel, EsmTokenizer

from pyeed.dbconnect import DatabaseConnector

logger = logging.getLogger(__name__)


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


def load_model_and_tokenizer(
    model_name: str,
) -> Tuple[
    Union[
        EsmModel,
        ESMC,
        ESM3InferenceClient,
        torch.nn.DataParallel,
        ESM3,
    ],  # Updated return type
    Union[EsmTokenizer, None],
    torch.device,
]:
    """
    Loads either an ESM++, ESM-3 (using ESMC or ESM3) or an ESM-2 (using Transformers) model,
    depending on the `model_name` provided. Uses multiple GPUs in parallel if available.

    Args:
        model_name (str): The model name or identifier.

    Returns:
        Tuple of (model, tokenizer, device)
    """
    token = get_hf_token()
    # Default device is the first CUDA device if available, else CPU.
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = None

    if "esmc" in model_name.lower():
        model: Any = ESMC.from_pretrained(model_name)
        model = model.to(device)
    elif "esm3-sm-open-v1" in model_name.lower():
        model: Any = ESM3.from_pretrained("esm3_sm_open_v1").to(device)
    else:
        full_model_name = (
            model_name
            if model_name.startswith("facebook/")
            else f"facebook/{model_name}"
        )
        model: Any = EsmModel.from_pretrained(full_model_name, use_auth_token=token)
        tokenizer = EsmTokenizer.from_pretrained(full_model_name, use_auth_token=token)
        model = model.to(device)

    # Check if multiple GPUs are available and wrap the model accordingly
    # if torch.cuda.device_count() > 1 and device.type == "cuda":
    #     logger.info(f"Using {torch.cuda.device_count()} GPUs for parallel inference.")
    #     model = torch.nn.DataParallel(model)

    return model, tokenizer, device


def get_batch_embeddings(
    batch_sequences: list[str],
    model: Union[
        EsmModel,
        ESMC,
        torch.nn.DataParallel,
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
    # base_model = model.module if isinstance(model, torch.nn.DataParallel) else model

    if isinstance(model, ESMC):
        # For ESMC models
        embedding_list = []
        with torch.no_grad():
            for sequence in batch_sequences:
                protein = ESMProtein(sequence=sequence)
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
                if pool_embeddings:
                    embeddings = embeddings.mean(axis=1)
                embedding_list.append(embeddings[0])
        return embedding_list
    elif isinstance(model, ESM3):
        # For ESM3 models
        embedding_list = []
        with torch.no_grad():
            for sequence in batch_sequences:
                protein = ESMProtein(sequence=sequence)
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
    sequence: str, model_name: str = "facebook/esm2_t33_650M_UR50D"
) -> NDArray[np.float64]:
    """
    Calculates an embedding for a single sequence.

    Args:
        sequence: Input protein sequence
        model_name: Name of the ESM model to use

    Returns:
        NDArray[np.float64]: Normalized embedding vector for the sequence
    """
    model, tokenizer, device = load_model_and_tokenizer(model_name)
    return get_single_embedding_last_hidden_state(sequence, model, tokenizer, device)


def calculate_single_sequence_embedding_all_layers(
    sequence: str, model_name: str = "facebook/esm2_t33_650M_UR50D"
) -> NDArray[np.float64]:
    """
    Calculates embeddings for a single sequence across all layers.

    Args:
        sequence: Input protein sequence
        model_name: Name of the ESM model to use

    Returns:
        NDArray[np.float64]: A numpy array containing layer embeddings for the sequence.
    """
    model, tokenizer, device = load_model_and_tokenizer(model_name)
    return get_single_embedding_all_layers(sequence, model, tokenizer, device)


def calculate_single_sequence_embedding_first_layer(
    sequence: str, model_name: str = "facebook/esm2_t33_650M_UR50D"
) -> NDArray[np.float64]:
    """
    Calculates an embedding for a single sequence using the first layer.
    """
    model, tokenizer, device = load_model_and_tokenizer(model_name)
    return get_single_embedding_first_layer(sequence, model, tokenizer, device)


def get_single_embedding_first_layer(
    sequence: str, model: Any, tokenizer: Any, device: torch.device
) -> NDArray[np.float64]:
    """
    Generates normalized embeddings for each token in the sequence across all layers.
    """
    embeddings_list = []

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
            if logits_output.hidden_states is None:
                raise ValueError(
                    "Model did not return hidden states. Check LogitsConfig settings."
                )
            embedding = (
                logits_output.hidden_states[0][0].to(torch.float32).cpu().numpy()
            )

        elif isinstance(model, ESM3):
            # ESM-3 logic
            from esm.sdk.api import ESMProtein, SamplingConfig

            protein = ESMProtein(sequence=sequence)
            protein_tensor = model.encode(protein)
            embedding = model.forward_and_sample(
                protein_tensor,
                SamplingConfig(return_per_residue_embeddings=True),
            )
            if embedding is None or embedding.per_residue_embedding is None:
                raise ValueError("Model did not return embeddings")
            embedding = embedding.per_residue_embedding.to(torch.float32).cpu().numpy()

        else:
            # ESM-2 logic
            inputs = tokenizer(sequence, return_tensors="pt").to(device)
            outputs = model(**inputs, output_hidden_states=True)
            # Get the first layer's hidden states for all residues (excluding special tokens)
            embedding = outputs.hidden_states[0][0, 1:-1, :].detach().cpu().numpy()

    # Ensure embedding is a numpy array and normalize it
    embedding = np.asarray(embedding, dtype=np.float64)
    embedding = embedding / np.linalg.norm(embedding, axis=1, keepdims=True)
    return embedding


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
                    sequence=True, return_embeddings=True, return_hidden_states=True
                ),
            )
            if logits_output.hidden_states is None:
                raise ValueError(
                    "Model did not return hidden states. Check LogitsConfig settings."
                )

            embedding = (
                logits_output.hidden_states[-1][0].to(torch.float32).cpu().numpy()
            )
        elif isinstance(model, ESM3):
            # ESM-3 logic
            from esm.sdk.api import ESMProtein, SamplingConfig

            protein = ESMProtein(sequence=sequence)
            sequence_encoding = model.encode(protein)

            embedding = model.forward_and_sample(
                sequence_encoding, SamplingConfig(return_per_residue_embeddings=True)
            )

            if embedding is None or embedding.per_residue_embedding is None:
                raise ValueError("Model did not return embeddings")
            embedding = embedding.per_residue_embedding.to(torch.float32).cpu().numpy()
        else:
            # ESM-2 logic
            inputs = tokenizer(sequence, return_tensors="pt").to(device)
            outputs = model(**inputs, output_hidden_states=True, return_dict=True)
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

        elif isinstance(model, ESM3):
            raise NotImplementedError("ESM3 is not supported for all layers")

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
