import logging
from typing import Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial as sp
from matplotlib.figure import Figure
from numpy.typing import NDArray
from scipy.spatial.distance import cosine

from pyeed.dbconnect import DatabaseConnector

logger = logging.getLogger(__name__)


class EmbeddingTool:
    def __init__(self) -> None:
        pass

    def get_embedding(
        self,
        sequence_id: str,
        db: DatabaseConnector,
    ) -> NDArray[np.float64]:
        """Get the embedding for a given sequence ID.

        Args:
            sequence_id (str): The sequence ID to retrieve the embedding for
            db (DatabaseConnector): The database connector object

        Returns:
            np.ndarray: The embedding vector

        Raises:
            ValueError: If no embedding is found for the given sequence ID
        """

        query = f"""
        MATCH (p:Protein {{accession_id: '{sequence_id}'}})
        RETURN p.embedding AS embedding
        """
        embedding_read = db.execute_read(query)
        if len(embedding_read) == 0:
            raise ValueError(f"No embedding found for sequence id: {sequence_id}")
        embedding = np.array(embedding_read[0]["embedding"])

        return embedding

    def find_closest_matches_simple(
        self,
        start_sequence_id: str,
        db: DatabaseConnector,
        metric: Literal["cosine", "euclidean", "manhattan"] = "cosine",
        n: int = 10,
    ) -> list[tuple[str, float]]:
        """Find the closest matches to the start_sequence_id in the embedding space.

        This function loads all embeddings into memory to calculate distances,
        which may not be optimal for large datasets.

        Args:
            start_sequence_id (str): The sequence ID to find matches for
            db (DatabaseConnector): The database connector object
            metric (str): The metric to use for distance calculation
            n (int): The number of closest matches to return

        Returns:
            list[tuple[str, float]]: List of tuples containing (sequence_id, distance)
        """

        # get the embedding for the start_sequence_id
        start_embedding = self.get_embedding(start_sequence_id, db)

        # get all the embeddings, where the embedding is not null
        query = """
        MATCH (p:Protein)
        WHERE p.embedding IS NOT NULL
        RETURN p.accession_id AS accession_id, p.embedding AS embedding
        """
        embeddings_read = db.execute_read(query)
        logger.info(f"Found {len(embeddings_read)} embeddings")

        # calculate the distances
        distances = []
        for embedding_read in embeddings_read:
            embedding = np.array(embedding_read["embedding"])
            if metric == "cosine":
                distance = cosine(start_embedding, embedding)
            elif metric == "euclidean":
                distance = np.linalg.norm(start_embedding - embedding)
            elif metric == "manhattan":
                distance = np.sum(np.abs(start_embedding - embedding))

            distances.append((embedding_read["accession_id"], distance))

        # sort the distances
        distances.sort(key=lambda x: x[1])

        return distances[:n]

    def calculate_2d_projection_tsne(
        self,
        db: DatabaseConnector,
        perplexity: int = 50,
        n_iter: int = 1000,
        ids_list: Optional[list[str]] = None,
        ids_list_labels: Optional[dict[str, str]] = None,
        random_state: int = 42,
    ) -> tuple[list[str], NDArray[np.float64], list[str], list[str]]:
        """Perform a 2D projection of the embeddings using t-SNE and prepare visualization data.

        Args:
            db (DatabaseConnector): The database connector object
            perplexity (int): The perplexity parameter for t-SNE
            n_iter (int): Number of iterations for t-SNE
            ids_list (Optional[list[str]]): List of sequence IDs to visualize
            ids_list_labels (Optional[dict[str, str]]): Dictionary mapping sequence IDs to labels
            random_state (int): Random seed for reproducibility

        Returns:
            tuple:
                - list[str]: List of protein IDs
                - np.ndarray: 2D projection coordinates
                - list[str]: Labels for each point
                - list[str]: Colors for each point
        """
        if ids_list is None:
            # get all the accession_ids
            query = """
            MATCH (p:Protein)
            WHERE p.accession_id IS NOT NULL
            RETURN p.accession_id AS protein_id
            """
            ids_list = [record["protein_id"] for record in db.execute_read(query)]

        # get the embeddings for the proteins based in the ids list
        query = """
        MATCH (p:Protein)
        WHERE p.accession_id IN $ids
        RETURN p.accession_id AS protein_id, p.embedding AS embedding
        """
        result = db.execute_read(query, {"ids": ids_list})

        # prepare data for visualization
        protein_ids, labels = [], []
        embeddings = np.array([record["embedding"] for record in result])

        for record in result:
            protein_ids.append(record["protein_id"])
            if ids_list_labels is not None:
                labels.append(ids_list_labels[record["protein_id"]])
            else:
                labels.append("None")

        # assign each label a color and create the color list with the corresponding colors
        # the default colo for 'None' is black
        colors = []
        color_label_dict: dict[str, str] = {}
        cycle_colors = [
            "blue",
            "red",
            "green",
            "orange",
            "purple",
            "brown",
            "pink",
            "gray",
            "olive",
            "cyan",
        ]

        for label in labels:
            if label not in color_label_dict.keys():
                color_label_dict[label] = cycle_colors[
                    len(color_label_dict) % len(cycle_colors)
                ]

        # assign the colors
        for label in labels:
            if label is None:
                colors.append("black")
            else:
                colors.append(color_label_dict[label])

        # perform t-SNE
        from sklearn.manifold import TSNE

        tsne = TSNE(
            n_components=2,
            perplexity=perplexity,
            max_iter=n_iter,
            random_state=random_state,
        )

        embeddings_2d = tsne.fit_transform(embeddings)

        return protein_ids, embeddings_2d, labels, colors

    def plot_matrix_comparison(
        self,
        distance_matrix_1: NDArray[np.float64],
        distance_matrix_2: NDArray[np.float64],
        protein_ids_1: list[str],
        protein_ids_2: list[str],
        label_1: str,
        label_2: str,
        number_of_points_goal: int,
    ) -> Figure:
        """Plot a comparison between two distance matrices.

        Args:
            distance_matrix_1 (NDArray[np.float64]): First distance matrix to compare
            distance_matrix_2 (NDArray[np.float64]): Second distance matrix to compare
            protein_ids_1 (list[str]): Protein IDs corresponding to first matrix
            protein_ids_2 (list[str]): Protein IDs corresponding to second matrix
            label_1 (str): Label for the first matrix
            label_2 (str): Label for the second matrix
            number_of_points_goal (int): Target number of points to plot

        Returns:
            plt.Figure: The generated plot
        """
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
                        plt.scatter(
                            distance_matrix_1[i, j],
                            distance_matrix_2[i, j],
                            c="gray",
                            alpha=0.1,
                            s=30,
                            edgecolor="k",
                        )
                    else:
                        # Convert the lists to np.array for valid elementwise comparison.
                        index_matrix_1_id = np.where(
                            np.array(protein_ids_1) == matrix_2_id
                        )[0][0]
                        index_matrix_2_id = np.where(
                            np.array(protein_ids_2) == matrix_1_id
                        )[0][0]

                        plt.scatter(
                            distance_matrix_1[i, index_matrix_1_id],
                            distance_matrix_2[index_matrix_2_id, j],
                            c="gray",
                            alpha=0.1,
                            s=30,
                            edgecolor="k",
                        )

                    point_counter += 1
                if point_counter > number_of_points_goal:
                    break

        plt.title(f"{label_1} vs {label_2}, points placed n = {point_counter}")
        plt.xlabel(label_1)
        plt.ylabel(label_2)
        plt.grid()
        plt.tight_layout()

        return fig

    def calculate_similarity(
        self,
        query_embed: NDArray[np.float64],
        target_embed: NDArray[np.float64],
        mode: Literal["cosine", "euclidean"] = "cosine",
    ) -> NDArray[np.float64]:
        """Calculate similarity between two protein sequence embeddings.

        Args:
            query_embed (NDArray[np.float64]): Embedding of the query sequence
            target_embed (NDArray[np.float64]): Embedding of the target sequence
            mode (Literal["cosine", "euclidean"]): Similarity metric to use

        Returns:
            np.ndarray: A 2D array containing similarity scores between all pairs
                of amino acids in the query and target sequences

        Raises:
            ValueError: If an invalid mode is specified
        """

        if mode == "cosine":
            return np.array(
                1 - sp.distance.cdist(query_embed, target_embed, metric="cosine")
            )
        elif mode == "euclidean":
            return np.array(
                1 / sp.distance.cdist(query_embed, target_embed, metric="euclidean")
            )
        else:
            raise ValueError(
                f"Invalid mode: {mode}, valid modes are 'cosine' or 'euclidean'"
            )

    def create_embedding_vector_index_neo4j(
        self,
        db: DatabaseConnector,
        index_name: str = "embedding_index",
        similarity_function: Literal["cosine", "euclidean"] = "cosine",
        dimensions: int = 1280,
        m: int = 16,
        ef_construction: int = 512,
    ) -> None:
        """
        Generates the native vector index for the embedding property in the Protein node. It uses the HNSW index type and the cosine distance metric.
        This is supported by Neo4j.

        The HNSW is a graph-based index that is optimized for high-dimensional data. It stands for Hierarchical Navigable Small World graphs.
        The the configuration the index is created. The parameters are:
        - index_type: "hnsw"
        - distance_metric: "cosine"
        - dimensions: 1280 (ESM2 model) and for ESM-C it is 960
        - name: "embedding_index"
        - m parameter: 16 it specifies the maximum number of outgoing edges per node
        - ef_construction parameter: 512 it specifies the size of the dynamic list for the nearest neighbors

        """

        query_create_index = f"""
            CREATE VECTOR INDEX {index_name}
            FOR (p:Protein) ON (p.embedding)
            OPTIONS {{
            indexProvider: 'vector-2.0',
            indexConfig: {{
                `vector.similarity_function`: '{similarity_function}',
                `vector.dimensions`: {dimensions},
                `vector.hnsw.m`: {m},
                `vector.hnsw.ef_construction`: {ef_construction},
                `vector.quantization.enabled`: False
            }}
            }};
        """
        db.execute_write(query_create_index)

    def find_nearest_neighbors_based_on_vector_index(
        self,
        query_id: str,
        db: DatabaseConnector,
        index_name: str = "embedding_index",
        number_of_neighbors: int = 50,
    ) -> list[tuple[str, float]]:
        """
        This function finds the nearest neighbors of a query protein based on the vector index.

        Args:
            db (DatabaseConnector): The database connector object
            index_name (str): The name of the vector index
            query_protein_id (str): The accession ID of the query protein
            number_of_neighbors (int): The number of nearest neighbors to find

        Returns:
            list[tuple[str, float]]: A list of tuples containing the accession ID and the similarity score
        """

        # Check if index is still being populated
        import time

        from rich.progress import Progress

        query_check_population = f"""
        SHOW INDEXES YIELD populationPercent, name 
        WHERE name = '{index_name}'
        RETURN populationPercent
        """

        # First check if index is already fully populated
        result = db.execute_read(query_check_population)
        if not result:
            raise ValueError(f"Index {index_name} not found")

        percent = result[0]["populationPercent"]
        if percent < 100:
            # If not fully populated, show progress bar
            with Progress() as progress:
                task = progress.add_task(
                    "[cyan]Waiting for index population...", total=100
                )
                progress.update(task, completed=percent)

                while True:
                    result = db.execute_read(query_check_population)
                    if not result:
                        raise ValueError(f"Index {index_name} not found")

                    percent = result[0]["populationPercent"]
                    if percent >= 100:
                        progress.update(task, completed=100)
                        break

                    progress.update(task, completed=percent)
                    time.sleep(0.01)
            logger.info(f"Index {index_name} is populated, finding nearest neighbors")

        query_find_nearest_neighbors = f"""
        MATCH (source:Protein {{accession_id: '{query_id}'}})
        WITH source.embedding AS embedding
        CALL db.index.vector.queryNodes('{index_name}', {number_of_neighbors}, embedding)
        YIELD node AS fprotein, score
        WHERE score > 0.95
        RETURN fprotein.accession_id, score
        """
        results = db.execute_read(query_find_nearest_neighbors)
        neighbors: list[tuple[str, float]] = [
            (str(record["fprotein.accession_id"]), float(record["score"]))
            for record in results
        ]
        return neighbors

    def drop_vector_index(
        self,
        db: DatabaseConnector,
        index_name: str = "embedding_index",
    ) -> None:
        """
        This function drops the vector index for the embedding property in the Protein node.

        Args:
            db (DatabaseConnector): The database connector object
            index_name (str): The name of the vector index
        """

        logger.info(f"Dropping vector index {index_name}")

        query_drop_index = f"DROP INDEX {index_name} IF EXISTS;"
        db.execute_write(query_drop_index)
