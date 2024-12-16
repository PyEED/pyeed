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

    def visualization_2d_projection_tsne(self, db: DatabaseConnector, perplexity = 50, n_iter = 1000, ids_list = None, ids_list_labels = None):
        """
        This function will perform a 2D projection of the embeddings using t-SNE and visualize it.
        
        Parameters:
        db (DatabaseConnector): The database connector object
        perplexity (int): The perplexity parameter for the t-SNE algorithm. Default is 30
        n_iter (int): The number of iterations for the t-SNE algorithm. Default is 1000
        ids_list (list): A list of sequence ids for which we want to visualize the embeddings. Default is None
        ids_list_labels (list): A list of labels for the sequence ids. Default is None
        """
        if ids_list is None:
            # get all the accession_ids
            query = """
            MATCH (p:Protein)
            WHERE p.accession_id IS NOT NULL
            RETURN p.accession_id AS protein_id
            """
            ids_list = [record['protein_id'] for record in db.execute_read(query)]

        
        # get the embeddings for the proteins based in the ids list
        query = """
        MATCH (p:Protein)
        WHERE p.accession_id IN $ids
        RETURN p.accession_id AS protein_id, p.embedding AS embedding
        """
        result = db.execute_read(query, {"ids": ids_list})
        print(f"Number of proteins fetched with embeddings: {len(result)}")
        print(f"Number of proteins in id list: {len(ids_list)}")

        # prepare data for visualization
        protein_ids, labels = [], []
        embeddings = []

        for record in result:
            protein_ids.append(record["protein_id"])
            embeddings.append(record["embedding"])
            # the label is either None or the label from the ids_list_labels
            if ids_list_labels is not None:
                labels.append(ids_list_labels[record["protein_id"]])
            else:
                labels.append('None')
            
        embeddings = np.array(embeddings)

        # assign each label a color and create the color list with the corresponding colors
        # the default colo for 'None' is black
        colors = []
        color_label_dict: dict[str, str] = {}
        import matplotlib.pyplot as plt
        cycle_colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']

        for label in labels:
            if label not in color_label_dict.keys():
                color_label_dict[label] = cycle_colors[len(color_label_dict) % len(cycle_colors)]
        
        # assign the colors
        for label in labels:
            if label == None:
                colors.append('black')
            else:
                colors.append(color_label_dict[label])
        
        # perform t-SNE
        from sklearn.manifold import TSNE
        tsne = TSNE(n_components = 2, perplexity = perplexity, max_iter = n_iter, random_state=42)

        embeddings_2d = tsne.fit_transform(embeddings)

        return protein_ids, embeddings_2d, labels, colors

    def visualization_2d_projection_umap(self, db: DatabaseConnector, n_neighbors = 15, min_dist = 0.1, metric = 'cosine', ids_list = None, ids_list_labels = None):
        """
        This function will perform a 2D projection of the embeddings using UMAP and visualize it.
        
        Parameters:
        db (DatabaseConnector): The database connector object
        n_neighbors (int): The number of neighbors parameter for the UMAP algorithm. Default is 15
        min_dist (float): The minimum distance parameter for the UMAP algorithm. Default is 0.1
        metric (str): The metric to use for the distance calculation. Default is 'cosine'
        ids_list (list): A list of sequence ids for which we want to visualize the embeddings. Default is None
        ids_list_labels (list): A list of labels for the sequence ids. Default is None
        """
        if ids_list is None:
            # get all the accession_ids
            query = """
            MATCH (p:Protein)
            WHERE p.accession_id IS NOT NULL
            RETURN p.accession_id AS protein_id
            """
            ids_list = [record['protein_id'] for record in db.execute_read(query)]

        
        # get the embeddings for the proteins based in the ids list
        query = """
        MATCH (p:Protein)
        WHERE p.accession_id IN $ids
        RETURN p.accession_id AS protein_id, p.embedding AS embedding
        """
        result = db.execute_read(query, {"ids": ids_list})
        print(f"Number of proteins fetched with embeddings: {len(result)}")
        print(f"Number of proteins in id list: {len(ids_list)}")

        # prepare data for visualization
        protein_ids, labels = [], []
        embeddings = []

        for record in result:
            protein_ids.append(record["protein_id"])
            print(type(embeddings))
            print(type(record["embedding"]))
            embeddings.append(record["embedding"])
            # the label is either None or the label from the ids_list_labels
            if ids_list_labels is not None:
                labels.append(ids_list_labels[record["protein_id"]])
            else:
                labels.append('None')
            
        embeddings = np.array(embeddings)

        # assign each label a color and create the color list with the corresponding colors
        # the default colo for 'None' is black
        colors = []
        color_label_dict: dict[str, str] = {}
        import matplotlib.pyplot as plt
        cycle_colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']


        for label in labels:
            if label not in color_label_dict.keys():
                color_label_dict[label] = cycle_colors[len(color_label_dict) % len(cycle_colors)]
            

        # assign the colors
        for label in labels:
            if label == None:
                colors.append('black')
            else:
                colors.append(color_label_dict[label])

        # perform UMAP
        import umap
        umap_model = umap.UMAP(n_neighbors = n_neighbors, min_dist = min_dist, metric = metric)
        embeddings_2d = umap_model.fit_transform(embeddings)

        return protein_ids, embeddings_2d, labels, colors


        

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


