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
    token = get_hf_token()
    tokenizer = None

    if "esmc" in model_name.lower():
        model = ESMC.from_pretrained(model_name)
        model = model.to(device)
    elif "esm3-sm-open-v1" in model_name.lower():
        model = ESM3.from_pretrained("esm3_sm_open_v1")
        model = model.to(device)
    elif "prot_t5" in model_name.lower() or "prott5" in model_name.lower():
        # ProtT5 models
        full_model_name = (
            model_name
            if model_name.startswith("Rostlab/")
            else f"Rostlab/{model_name}"
        )
        model = T5Model.from_pretrained(full_model_name, use_auth_token=token)
        tokenizer = T5Tokenizer.from_pretrained(full_model_name, use_auth_token=token, do_lower_case=False)
        model = model.to(device)
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
    elif isinstance(base_model, T5Model):
        # For ProtT5 models
        assert tokenizer_or_alphabet is not None, "Tokenizer required for ProtT5 models"
        assert isinstance(tokenizer_or_alphabet, T5Tokenizer), "T5Tokenizer required for ProtT5 models"
        
        # Preprocess sequences for ProtT5
        processed_sequences = [preprocess_sequence_for_prott5(seq) for seq in batch_sequences]
        
        inputs = tokenizer_or_alphabet.batch_encode_plus(
            processed_sequences, 
            add_special_tokens=True, 
            padding="longest",
            return_tensors="pt"
        )
        
        # Move inputs to device
        input_ids = inputs['input_ids'].to(device)
        attention_mask = inputs['attention_mask'].to(device)
        
        with torch.no_grad():
            # For ProtT5, use encoder embeddings for feature extraction
            # Create dummy decoder inputs (just the pad token)
            batch_size = input_ids.shape[0]
            decoder_input_ids = torch.full(
                (batch_size, 1), 
                tokenizer_or_alphabet.pad_token_id or 0, 
                dtype=torch.long,
                device=device
            )
            
            outputs = base_model(input_ids=input_ids, 
                          attention_mask=attention_mask,
                          decoder_input_ids=decoder_input_ids)
            
            # Get encoder last hidden state (encoder embeddings)
            hidden_states = outputs.encoder_last_hidden_state.cpu().numpy()

        if pool_embeddings:
            # Mean pooling across sequence length, excluding padding tokens
            embedding_list = []
            for i, hidden_state in enumerate(hidden_states):
                # Get actual sequence length (excluding padding)
                attention_mask_np = attention_mask[i].cpu().numpy()
                seq_len = attention_mask_np.sum()
                # Pool only over actual sequence tokens
                pooled_embedding = hidden_state[:seq_len].mean(axis=0)
                embedding_list.append(pooled_embedding)
            return embedding_list
        return list(hidden_states)
    else:
        # ESM-2 logic
        assert tokenizer_or_alphabet is not None, "Tokenizer required for ESM-2 models"
        assert isinstance(tokenizer_or_alphabet, EsmTokenizer), "EsmTokenizer required for ESM-2 models"
        inputs = tokenizer_or_alphabet(
            batch_sequences, padding=True, truncation=True, return_tensors="pt"
        ).to(device)
        with torch.no_grad():
            outputs = base_model(**inputs, output_hidden_states=True)

        # Get last hidden state for each sequence
        hidden_states = outputs.last_hidden_state.cpu().numpy()

        if pool_embeddings:
            # Mean pooling across sequence length
            return [embedding.mean(axis=0) for embedding in hidden_states]
        return list(hidden_states)


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
        elif isinstance(model, T5Model):
            # ProtT5 logic
            processed_sequence = preprocess_sequence_for_prott5(sequence)
            inputs = tokenizer.encode_plus(
                processed_sequence,
                add_special_tokens=True,
                return_tensors="pt"
            )
            
            input_ids = inputs['input_ids'].to(device)
            attention_mask = inputs['attention_mask'].to(device)
            
            # Create dummy decoder inputs
            decoder_input_ids = torch.full(
                (1, 1), 
                tokenizer.pad_token_id or 0, 
                dtype=torch.long,
                device=device
            )
            
            outputs = model(input_ids=input_ids, 
                          attention_mask=attention_mask,
                          decoder_input_ids=decoder_input_ids)
            
            # Get encoder last hidden state including special tokens
            embedding = outputs.encoder_last_hidden_state[0].detach().cpu().numpy()
        else:
            # ESM-2 logic
            inputs = tokenizer(sequence, return_tensors="pt").to(device)
            outputs = model(**inputs)
            embedding = outputs.last_hidden_state[0, 1:-1, :].detach().cpu().numpy()

    return embedding  # type: ignore


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

        elif isinstance(model, T5Model):
            # For ProtT5: Get encoder hidden states
            processed_sequence = preprocess_sequence_for_prott5(sequence)
            inputs = tokenizer.encode_plus(
                processed_sequence,
                add_special_tokens=True,
                return_tensors="pt"
            )
            
            input_ids = inputs['input_ids'].to(device)
            attention_mask = inputs['attention_mask'].to(device)
            
            # Create dummy decoder inputs
            decoder_input_ids = torch.full(
                (1, 1), 
                tokenizer.pad_token_id or 0, 
                dtype=torch.long,
                device=device
            )
            
            outputs = model(input_ids=input_ids, 
                          attention_mask=attention_mask,
                          decoder_input_ids=decoder_input_ids,
                          output_hidden_states=True)
            
            # Get all encoder hidden states
            encoder_hidden_states = outputs.encoder_hidden_states
            for layer_tensor in encoder_hidden_states:
                # Remove batch dimension but keep special tokens
                emb = layer_tensor[0].detach().cpu().numpy()
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

def calculate_single_sequence_embedding_first_layer(
    sequence: str, model_name: str = "facebook/esm2_t33_650M_UR50D", device: torch.device = torch.device("cuda:0"),
) -> NDArray[np.float64]:
    """
    Calculates an embedding for a single sequence using the first layer.
    """
    model, tokenizer, device = load_model_and_tokenizer(model_name, device)
    return get_single_embedding_first_layer(sequence, model, tokenizer, device)


def get_single_embedding_first_layer(
    sequence: str, model: Any, tokenizer: Any, device: torch.device
) -> NDArray[np.float64]:
    """
    Generates normalized embeddings for each token in the sequence using the first layer.
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

        elif isinstance(model, T5Model):
            # ProtT5 logic - get first layer embedding
            processed_sequence = preprocess_sequence_for_prott5(sequence)
            inputs = tokenizer.encode_plus(
                processed_sequence,
                add_special_tokens=True,
                return_tensors="pt"
            )
            
            input_ids = inputs['input_ids'].to(device)
            attention_mask = inputs['attention_mask'].to(device)
            
            # Create dummy decoder inputs
            decoder_input_ids = torch.full(
                (1, 1), 
                tokenizer.pad_token_id or 0, 
                dtype=torch.long,
                device=device
            )
            
            outputs = model(input_ids=input_ids, 
                          attention_mask=attention_mask,
                          decoder_input_ids=decoder_input_ids,
                          output_hidden_states=True)
            
            # Get first encoder hidden state including special tokens
            embedding = outputs.encoder_hidden_states[0][0].detach().cpu().numpy()

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
