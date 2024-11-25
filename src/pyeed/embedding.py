import gc
import logging

import torch
from transformers import EsmModel, EsmTokenizer

from pyeed.dbconnect import DatabaseConnector

logger = logging.getLogger(__name__)


def load_model_and_tokenizer(model_name: str):
    """
    Loads the ESM2 model and tokenizer and sets the appropriate device.
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


def get_batch_embeddings(batch_sequences, model, tokenizer, device):
    """
    Generates embeddings for a batch of sequences.
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
    embeddings_batch,
):
    """
    Updates the embeddings for a batch of proteins in the database.
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


def free_memory():
    """
    Frees up memory by invoking garbage collection and clearing GPU caches.
    """
    gc.collect()
    if torch.backends.mps.is_available():
        torch.mps.empty_cache()
    elif torch.cuda.is_available():
        torch.cuda.empty_cache()
