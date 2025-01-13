# this is basically a tool to analyze the embedding results
# it assumes that the embedding is already done and the results are saved in the neo4j database
# it will load the embedding results and perform some analysis on it
# a lot of code was provided by Tim Panzer in his bachelor thesis results
import torch
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial as sp

from pyeed.main import Pyeed
from pyeed.dbconnect import DatabaseConnector
from pyeed.embedding import load_model_and_tokenizer


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
        if len(embedding_read) == 0:
            raise ValueError(f"No embedding found for sequence id: {sequence_id}")
        embedding = np.array(embedding_read[0]['embedding'])

        return embedding

    def _get_single_embedding_last_hidden_state(self, sequence, model, tokenizer, device):
        """
        Generates embeddings for a single sequence.
        And return the last hidden state. As a numpy array. Not the mean of the last hidden state.
        Allows analysis of the single token in the sequences.
        """

        with torch.no_grad():
            inputs = tokenizer(sequence, return_tensors="pt").to(device)
            outputs = model(**inputs)
            embedding = outputs.last_hidden_state[0, 1:-1, :].detach().cpu().numpy()

        # normalize the embedding
        embedding = embedding / np.linalg.norm(embedding, axis=1, keepdims=True)

        return embedding

    def calculate_single_sequence_embedding(self, sequence: str, model_name: str = "facebook/esm2_t33_650M_UR50D"):
        """
        Calculates an embedding for a single sequence.
        """
        model, tokenizer, device = load_model_and_tokenizer(model_name)
        return self._get_single_embedding_last_hidden_state(sequence, model, tokenizer, device)

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

    def plot_matrix_comparison(self, distance_matrix_1, distance_matrix_2, protein_ids_1, protein_ids_2, label_1, label_2, number_of_points_goal):
        # general plot function
        # we want to plot the two distance matrices against each other
        fig = plt.figure(figsize=(15, 10))

        total_number_of_points = distance_matrix_1.shape[0] * distance_matrix_1.shape[1]
        random_chance_probability = number_of_points_goal / total_number_of_points
        point_counter = 0

        for i in range(len(protein_ids_1)):
            for j in range(len(protein_ids_2)):
                if np.random.rand() < random_chance_probability:

                    matrix_1_id = protein_ids_1[i]
                    matrix_2_id = protein_ids_2[j]

                    if matrix_1_id == matrix_2_id:
                        plt.scatter(distance_matrix_1[i, j], distance_matrix_2[i, j], c='gray', alpha=0.1, s=30, edgecolor="k")
                    else:
                        # we need to find the index of the matrix_2_id in the protein_ids_1 array from numpy
                        # the i is the index in the protein_ids_1 array this means we are looking of the index to j in the protein_ids_1 array
                        index_matrix_1_id = np.where(protein_ids_1 == matrix_2_id)[0][0]
                        # now the index of the matrix_1_id is the index of the matrix_2_id in the protein_ids_2 array
                        index_matrix_2_id = np.where(protein_ids_2 == matrix_1_id)[0][0]

                        plt.scatter(distance_matrix_1[i, index_matrix_1_id], distance_matrix_2[index_matrix_2_id, j], c='gray', alpha=0.1, s=30, edgecolor="k")
                    
                    point_counter += 1
                if point_counter > number_of_points_goal:
                    break
        
        plt.title(f"{label_1} vs {label_2}, points placed n = {point_counter}")
        plt.xlabel(label_1)
        plt.ylabel(label_2)
        plt.grid()
        plt.tight_layout()
        plt.show()

    def calculate_similarity(
        self, 
        query_embed: np.ndarray,
        target_embed: np.ndarray,
        mode: Literal["cosine", "euclidean"] = "cosine",
    ) -> np.ndarray:
        """Calculate cosine or euclidean similarity between two protein sequences.

        Args:
            query_embed (np.ndarray): Embedding of the query sequence
            target_embed (np.ndarray): Embedding of the target sequence
            mode (str): The mode of similarity to calculate. Can be "cosine" or "euclidean".
        Returns:
            numpy.ndarray: A 2D numpy array containing cosine similarity scores
                between all pairs of amino acids in seq1 and seq2.
        """

        if mode == "cosine":
            return 1 - sp.distance.cdist(query_embed, target_embed, metric="cosine")
        elif mode == "euclidean":
            return 1 / sp.distance.cdist(query_embed, target_embed, metric="euclidean")
        else:
            raise ValueError(f"Invalid mode: {mode}")
        

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


