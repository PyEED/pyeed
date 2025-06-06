"""
Database operations for protein embeddings.

Handles storing and updating protein embeddings in the database.
"""

from typing import List

import numpy as np
from numpy.typing import NDArray

from pyeed.dbconnect import DatabaseConnector


def update_protein_embeddings_in_db(
    db: DatabaseConnector,
    accessions: List[str],
    embeddings_batch: List[NDArray[np.float64]],
) -> None:
    """
    Updates the embeddings for a batch of proteins in the database.

    Args:
        db (DatabaseConnector): The database connector.
        accessions (List[str]): The accessions of the proteins to update.
        embeddings_batch (List[NDArray[np.float64]]): The embeddings to update.
    """
    # Prepare the data for batch update
    updates = []
    for acc, emb in zip(accessions, embeddings_batch):
        # Flatten the embedding array and convert to list
        flat_embedding = emb.flatten().tolist()
        updates.append({"accession": acc, "embedding": flat_embedding})

    # Cypher query for batch update
    query = """
    UNWIND $updates AS update
    MATCH (p:Protein {accession_id: update.accession})
    SET p.embedding = update.embedding
    """

    # Execute the update query with parameters
    db.execute_write(query, {"updates": updates})
