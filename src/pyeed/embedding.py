import gc
import os
from typing import Any, Tuple, Union

import numpy as np
import torch
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig
from huggingface_hub import HfFolder, login
from numpy.typing import NDArray
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


def load_model_and_tokenizer(
    model_name: str,
) -> Tuple[
    Union[EsmModel, ESMC],  # Changed from ESM3InferenceClient to ESMC
    Union[EsmTokenizer, None],
    torch.device,
]:
    """
    Loads either an ESM-3 (using ESMC) or an ESM-2 (using Transformers) model,
    depending on the `model_name` provided.

    Args:
        model_name (str): The model name or identifier (e.g., 'esmc' or 'esm2_t12_35M_UR50D').

    Returns:
        Tuple of (model, tokenizer, device)
    """
    # Get token only when loading model
    token = get_hf_token()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Check if this is an ESM-3 variant
    if "esmc" in model_name.lower():
        # Using ESMC from_pretrained
        model = ESMC.from_pretrained(model_name)
        model = model.to(device)
        return model, None, device
    else:
        # Otherwise, assume it's an ESM-2 model on Hugging Face
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
    model: Union[EsmModel, ESMC],  # Updated type hint
    tokenizer_or_alphabet: Union[EsmTokenizer, None],
    device: torch.device,
    pool_embeddings: bool = True,
) -> list[NDArray[np.float64]]:
    """
    Generates mean-pooled embeddings for a batch of sequences.

    Args:
        batch_sequences (list[str]): List of sequence strings to be embedded.
        model (Union[EsmModel, ESMC]): Loaded model (ESM-2 or ESM-3).
        tokenizer_or_alphabet (Union[EsmTokenizer, None]): Tokenizer if ESM-2, None if ESM-3.
        device (torch.device): Device on which to run inference (CPU or GPU).
        pool_embeddings (bool): Whether to pool embeddings across sequence length.

    Returns:
        list[NDArray[np.float64]]: A list of embeddings as NumPy arrays.
    """
    if isinstance(model, ESMC):
        with torch.no_grad():
            embedding_list = []
            for sequence in batch_sequences:
                # Process each sequence individually
                protein = ESMProtein(sequence=sequence)
                protein_tensor = model.encode(protein)
                logits_output = model.logits(
                    protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
                )
                # Convert embeddings to numpy array
                embeddings = logits_output.embeddings.cpu().numpy()
                if pool_embeddings:
                    embeddings = embeddings.mean(axis=1)
                embedding_list.append(embeddings[0])

        return embedding_list

    else:
        # ESM-2 logic
        assert tokenizer_or_alphabet is not None, "Tokenizer required for ESM-2 models"
        inputs = tokenizer_or_alphabet(
            batch_sequences, padding=True, truncation=True, return_tensors="pt"
        ).to(device)
        with torch.no_grad():
            outputs = model(**inputs)
        embeddings = outputs.last_hidden_state.cpu().numpy()
        if pool_embeddings:
            return [embedding.mean(axis=0) for embedding in embeddings]
        return list(embeddings)


def calculate_single_sequence_embedding_last_hidden_state(
    sequence: str, model_name: str = "facebook/esm2_t33_650M_UR50D"
):
    """
    Calculates an embedding for a single sequence.
    """
    model, tokenizer, device = load_model_and_tokenizer(model_name)
    return get_single_embedding_last_hidden_state(sequence, model, tokenizer, device)


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
                protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
            )
            embedding = logits_output.embeddings[0].cpu().numpy()
        else:
            # ESM-2 logic
            inputs = tokenizer(sequence, return_tensors="pt").to(device)
            outputs = model(**inputs)
            embedding = outputs.last_hidden_state[0, 1:-1, :].detach().cpu().numpy()

    # normalize the embedding
    embedding = embedding / np.linalg.norm(embedding, axis=1, keepdims=True)

    return embedding  # type: ignore


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
