# this is basically a tool to analyze the embedding results
# it assumes that the embedding is already done and the results are saved in the neo4j database
# it will load the embedding results and perform some analysis on it
# a lot of code was provided by Tim Panzer in his bachelor thesis results
import numpy as np
from scipy.spatial.distance import cosine

from pyeed.main import Pyeed
from pyeed.dbconnect import DatabaseConnector



class EmbeddingTool:
    def __init__(self):
        None

    def get_embedding(self, sequence_id: str, db: DatabaseConnector):
        """
        This function will get the embedding for a given sequence id
        
        Parameters:
        sequence_id (str): The sequence id for which we want to get the embedding
        db (DatabaseConnector): The database connector object

        Returns:
        numpy.ndarray: The embedding
        """

        query = f"""
        MATCH (p:Protein {{accession_id: '{sequence_id}'}})
        RETURN p.embedding AS embedding
        """
        embedding_read = db.execute_read(query)
        embedding = np.array(embedding_read[0]['embedding'])

        return embedding

    def find_closest_matches_simple(self, start_sequence_id: str, db: DatabaseConnector, metric = 'cosine', n = 10):
        """
        This function will find the closest matches to the start_sequence_id in the embedding space.
        The simple in the name indicates that this function will load all the embeddings into memory and calculate the distances.
        Hence being simple but maybe quite bad in terms of performance. ;)
        
        Parameters:
        start_sequence_id (str): The sequence id for which we want to find the closest matches
        db (DatabaseConnector): The database connector object
        metric (str): The metric to use for the distance calculation. Default is 'cosine'
        n (int): The number of closest matches to return. Default is 10

        Returns:
        list: A list of tuples containing the sequence id and the distance
        """

        # get the embedding for the start_sequence_id
        start_embedding = self.get_embedding(start_sequence_id, db)

        # get all the embeddings
        query = """
        MATCH (p:Protein)
        RETURN p.accession_id AS accession_id, p.embedding AS embedding
        """
        embeddings_read = db.execute_read(query)

        # calculate the distances
        distances = []
        for embedding_read in embeddings_read:
            embedding = np.array(embedding_read['embedding'])
            if metric == 'cosine':
                distance = cosine(start_embedding, embedding)
            elif metric == 'euclidean':
                distance = np.linalg.norm(start_embedding - embedding)
            elif metric == 'manhattan':
                distance = np.sum(np.abs(start_embedding - embedding))
            
            distances.append((embedding_read['accession_id'], distance))

        # sort the distances
        distances.sort(key = lambda x: x[1])

        return distances[:n]



if __name__ == "__main__":

    uri = "bolt://127.0.0.1:8123"
    user = "neo4j"
    password = "niklasniklaspwtem"
    
    eedb = Pyeed(uri, user=user, password=password)

    et = EmbeddingTool()

    embedding_single = et.get_embedding("AQT03459.1", eedb.db)
    print(f"Embedding for AQT03459.1: {embedding_single}")

    closest_matches = et.find_closest_matches_simple("AQT03459.1", eedb.db, metric = 'cosine', n = 10)
    print(f"Closest matches for AQT03459.1: {closest_matches}")







