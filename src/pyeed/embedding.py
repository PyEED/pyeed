import gc

import numpy as np
import torch
from numpy.typing import NDArray
from transformers import EsmModel, EsmTokenizer

from pyeed.dbconnect import DatabaseConnector


def load_model_and_tokenizer(
    model_name: str,
) -> tuple[EsmModel, EsmTokenizer, torch.device]:
    """
    Loads the ESM2 model and tokenizer and sets the appropriate device.

    Args:
        model_name (str): The name of the ESM2 model to load.

    Returns:
        tuple[EsmModel, EsmTokenizer, torch.device]: The ESM2 model, tokenizer, and device.
    """
    # Load the model and tokenizer
    model = EsmModel.from_pretrained(model_name)
    tokenizer = EsmTokenizer.from_pretrained(model_name)

    # Determine the device
    if torch.backends.mps.is_available():
        device = torch.device("mps")
    elif torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")

    model = model.to(device)
    model.eval()

    return model, tokenizer, device


def get_batch_embeddings(
    batch_sequences: list[str],
    model: EsmModel,
    tokenizer: EsmTokenizer,
    device: torch.device,
) -> list[NDArray[np.float64]]:
    """
    Generates embeddings for a batch of sequences.

    Args:
        batch_sequences (list[str]): The sequences to embed.
        model (EsmModel): The ESM2 model.
        tokenizer (EsmTokenizer): The ESM2 tokenizer.
        device (torch.device): The device to use for the embeddings.

    Returns:
        list[NDArray[np.float64]]: The embeddings for the sequences.
    """
    max_sequence_length = 1024
    with torch.no_grad():
        # Tokenize the input sequences
        inputs = tokenizer(
            batch_sequences,
            padding=True,
            truncation=True,
            return_tensors="pt",
            max_length=max_sequence_length,
        ).to(device)

        # Get model outputs
        outputs = model(**inputs)
        embeddings = outputs.last_hidden_state

        embedding_list = []
        # Process each sequence in the batch
        for j in range(len(batch_sequences)):
            valid_token_mask = inputs["attention_mask"][j].bool()
            seq_embeddings = embeddings[j][valid_token_mask].mean(dim=0).cpu()
            embedding_list.append(seq_embeddings.numpy())

    return embedding_list


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


def free_memory() -> None:
    """
    Frees up memory by invoking garbage collection and clearing GPU caches.
    """
    gc.collect()
    if torch.backends.mps.is_available():
        torch.mps.empty_cache()
    elif torch.cuda.is_available():
        torch.cuda.empty_cache()


if __name__ == "__main__":
    model_name = "esmc_300m"

    model, tokenizer, device = load_model_and_tokenizer(model_name)

    print(model)
    print(tokenizer)
    print(device)
