from typing import Any

import networkx as nx
from loguru import logger
from pyeed.dbconnect import DatabaseConnector


class NetworkAnalysis:
    """
    Class for network analysis and graph creation based on the Neo4j model.
    """

    def __init__(self, db: DatabaseConnector):
        """
        Initializes the NetworkAnalysis class.

        Args:
            db (DatabaseConnector): The database connector object.
        """
        self.db: DatabaseConnector = db
        self.graph: nx.Graph = nx.Graph()

    def _build_node_filter(
        self, node_types: list[str], ids: list[str]
    ) -> tuple[str, dict[str, Any]]:
        """Build the node filter query and parameters."""
        node_filter = ""
        params = {}

        if node_types:
            node_filter += "WHERE labels(n)[0] IN $node_types "
            params["node_types"] = node_types

        if ids:
            if node_filter:
                node_filter += "AND n.accession_id IN $accession_ids "
            else:
                node_filter += "WHERE n.accession_id IN $accession_ids "
            params["accession_ids"] = ids

        return node_filter, params

    def _fetch_nodes(
        self, node_filter: str, params: dict[str, Any]
    ) -> list[dict[str, Any]]:
        """Fetch nodes from the database."""
        query = f"""
        MATCH (n)
        {node_filter}
        RETURN ID(n) as id, labels(n) as labels, properties(n) as properties
        """
        logger.debug(f"Executing query: {query}")
        return self.db.execute_read(query, params)

    def _fetch_relationships(self, relationships: list[str]) -> list[dict[str, Any]]:
        """Fetch relationships from the database."""
        rel_filter = "WHERE type(r) IN $relationships " if relationships else ""
        query = f"""
        MATCH (n)-[r]->(m)
        {rel_filter}
        RETURN ID(n) as source, ID(m) as target, type(r) as type, properties(r) as properties
        """
        logger.debug(f"Executing query: {query}")
        return self.db.execute_read(query, {"relationships": relationships})

    def create_graph(
        self,
        node_types: list[str] = [],
        relationships: list[str] = [],
        ids: list[str] = [],
    ) -> nx.Graph:
        """Creates a graph using data from the Neo4j database with specified filters."""
        logger.info(
            f"Creating graph with node types: {node_types}, relationships: {relationships}, ids: {ids}"
        )

        # Fetch data
        node_filter, params = self._build_node_filter(node_types, ids)
        nodes = self._fetch_nodes(node_filter, params)
        relationships_result = self._fetch_relationships(relationships)

        # Build graph
        self.graph = nx.Graph()

        # Add nodes
        for node in nodes:
            self.graph.add_node(
                node["id"],
                labels=node["labels"],
                properties=node["properties"],
            )

        # Add relationships
        for rel in relationships_result:
            if rel["source"] in self.graph and rel["target"] in self.graph:
                self.graph.add_edge(
                    rel["source"],
                    rel["target"],
                    type=rel["type"],
                    properties=rel["properties"],
                )

        logger.debug(
            f"Created graph with {len(nodes)} nodes and {len(relationships_result)} relationships"
        )
        return self.graph

    def compute_degree_centrality(self) -> dict[Any, float]:
        """
        Computes the degree centrality of the graph.

        Returns:
            dict[Any, float]: A dictionary mapping each node to its degree centrality.
        """
        return dict(nx.degree_centrality(self.graph))

    def detect_communities(self) -> list[list[Any]]:
        """
        Detects communities in the graph using the greedy modularity algorithm.

        Returns:
            list[list[Any]]: A list of lists, each containing the nodes in a community.
        """
        from networkx.algorithms.community import greedy_modularity_communities

        # Convert frozenset to list explicitly
        return [list(c) for c in greedy_modularity_communities(self.graph)]

    def find_isolated_nodes(self, graph: nx.Graph) -> list[Any]:
        """
        Finds isolated nodes in the graph.

        Args:
            graph (networkx.Graph): The graph to analyze.

        Returns:
            list[Any]: List of nodes that have no connections (degree = 0).
        """
        degrees = dict(graph.degree())
        return [node for node in graph.nodes() if degrees[node] == 0]

    def find_self_referential_nodes(
        self, relationship_type: str | None, graph: nx.Graph
    ) -> list[tuple[Any, Any, dict[str, Any]]]:
        """
        Finds all edges in the graph that are self-referential (loops) of a specified relationship type.

        Args:
            relationship_type (str | None): The relationship type to consider. Can be None to match all types.
            graph (networkx.Graph): The graph to analyze.

        Returns:
            list[tuple]: A list of tuples representing the self-referential edges,
                        where each tuple is (node, node, edge_data).
        """
        removed_edges = []

        for node in list(graph.nodes):  # Use list to avoid mutation issues
            # Get all edges connected to the node
            edges = [
                (u, v, d)
                for u, v, d in graph.edges(node, data=True)
                if d.get("type") == relationship_type
            ]

            # Separate self-referential edges and others
            self_edges = [(u, v, d) for u, v, d in edges if u == v]
            other_edges = [(u, v, d) for u, v, d in edges if u != v]

            # If there are other edges, remove self-referential ones
            if other_edges and self_edges:
                for u, v, d in self_edges:
                    removed_edges.append((u, v, d))

        return removed_edges

    def calculate_positions_2d(
        self,
        attribute: str | None = None,
        scale: int = 1,
        threshold: float | None = None,
        mode: str = "HIDE_UNDER_THRESHOLD",
        type_relationship: str | None = None,
    ) -> tuple[nx.Graph, dict[Any, tuple[float, float]]]:
        """
        Calculates 2D positions for graph visualization with optional edge filtering.

        Args:
            attribute (str | None): Edge attribute to use for weighting. Can be None for unweighted layout.
            scale (int): Scale factor for node distances in the visualization. Defaults to 1.
            threshold (float | None): Threshold value for edge filtering.
            mode (str): Filtering mode - either "HIDE_UNDER_THRESHOLD" or "HIDE_OVER_THRESHOLD".
                       Defaults to "HIDE_UNDER_THRESHOLD".
            type_relationship (str | None): Specific relationship type to filter on. Can be None to match all types.

        Returns:
            tuple[nx.Graph, dict[Any, tuple[float, float]]]: A tuple containing:
                - The filtered graph
                - Dictionary mapping node IDs to their 2D positions as (x, y) coordinates
        """

        # Filter edges based on the threshold
        edges = self.graph.edges(data=True)
        edges_filtered = []

        if threshold is not None and attribute is not None:
            for edge in edges:
                if mode == "HIDE_UNDER_THRESHOLD":
                    if type_relationship is not None:
                        if edge[2]["type"] == type_relationship:
                            if edge[2]["properties"][attribute] >= threshold:
                                edges_filtered.append(edge)
                    else:
                        if edge[2]["properties"][attribute] >= threshold:
                            edges_filtered.append(edge)
                elif mode == "HIDE_OVER_THRESHOLD":
                    if type_relationship is not None:
                        if edge[2]["type"] == type_relationship:
                            if edge[2]["properties"][attribute] <= threshold:
                                edges_filtered.append(edge)
                    else:
                        if edge[2]["properties"][attribute] <= threshold:
                            edges_filtered.append(edge)

        logger.debug(f"Number of edges filtered: {len(edges) - len(edges_filtered)}")

        # Extract filtered graph
        filtered_graph = nx.Graph()
        filtered_graph.add_nodes_from(self.graph.nodes(data=True))
        filtered_graph.add_edges_from(edges_filtered)

        # Find self-referential nodes
        self_referential_edges = self.find_self_referential_nodes(
            type_relationship, filtered_graph
        )
        logger.info(f"Number of self-referential edges: {len(self_referential_edges)}")
        filtered_graph.remove_edges_from(self_referential_edges)

        # Find isolated nodes
        isolated_nodes = self.find_isolated_nodes(filtered_graph)
        filtered_graph.remove_nodes_from(isolated_nodes)

        # Use spring layout for force-directed graph
        weight_attr = attribute if attribute is not None else None
        pos = nx.spring_layout(filtered_graph, weight=weight_attr, scale=scale)
        # Convert pos to the correct return type
        positions: dict[Any, tuple[float, float]] = {
            node: (coord[0], coord[1]) for node, coord in pos.items()
        }

        return filtered_graph, positions

    def compute_clustering_coefficients(self) -> tuple[dict[Any, float], float]:
        """
        Computes the clustering coefficient for all nodes and the overall graph.

        Returns:
            tuple[dict[Any, float], float]: A tuple containing:
                - Dictionary mapping nodes to their clustering coefficients
                - Overall graph clustering coefficient
        """
        clustering = nx.clustering(self.graph)
        node_coefficients: dict[Any, float] = {
            node: float(coeff) for node, coeff in clustering.items()
        }
        avg_clustering = float(nx.average_clustering(self.graph))
        return node_coefficients, avg_clustering

    def analyze_connectivity(self) -> dict[str, Any]:
        """
        Analyzes the connectivity of the graph.

        Returns:
            dict[str, Any]: Information about whether the graph is connected and its connected components.
        """
        is_connected = nx.is_connected(self.graph)
        components = list(nx.connected_components(self.graph))
        return {"is_connected": is_connected, "components": components}

    def analyze_node(self, node_id: Any) -> dict[str, Any]:
        """
        Analyzes the properties and relationships of a specific node.

        Args:
            node_id (Any): The ID of the node to analyze.

        Returns:
            dict[str, Any]: Information about the node including:
                - id: The node ID
                - labels: Node labels from the graph
                - properties: Node properties
                - degree: Number of connections
                - neighbors: List of neighboring nodes
        """
        if node_id not in self.graph:
            raise ValueError(f"Node {node_id} does not exist in the graph.")

        node_data = self.graph.nodes[node_id]
        neighbors = list(self.graph.neighbors(node_id))
        degree = self.graph.degree(node_id)

        return {
            "id": node_id,
            "labels": node_data.get("labels"),
            "properties": node_data.get("properties"),
            "degree": degree,
            "neighbors": neighbors,
        }


if __name__ == "__main__":
    import logging

    from pyeed import Pyeed
    from pyeed.analysis.sequence_alignment import PairwiseAligner

    # Set up logging
    logging.basicConfig(
        level=logging.ERROR, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Database connection
    uri = "bolt://127.0.0.1:7688"
    user = "neo4j"
    password = "12345678"

    eedb = Pyeed(uri, user=user, password=password)

    # Example IDs for testing
    test_ids = [
        "2FX5_A",
        "ACC95208.1",
        "ACY95991.1",
        "ACY96861.1",
        "ADH43200.1",
        "ADK73612.1",
        "ADM47605.1",
        "ADV92525.1",
        "ADV92526.1",
        "ADV92528.1",
        "ADV92571.1",
        "AET05798.1",
        # ... add more IDs as needed
    ]

    # Wipe database and remove constraints
    # eedb.db.wipe_database()

    # Fetch sequences
    eedb.fetch_from_primary_db(test_ids, db="ncbi_protein")

    # Perform sequence alignment
    pa = PairwiseAligner()
    pa.align_multipairwise(ids=test_ids, db=eedb.db)

    # Create and analyze network
    na = NetworkAnalysis(db=eedb.db)
    na.create_graph(ids=test_ids, node_types=["Protein"])

    # Calculate positions with filtering
    attribute = "similarity"
    scale = 1
    threshold = 0.15
    mode = "HIDE_UNDER_THRESHOLD"
    type_relationship = "PAIRWISE_ALIGNED"

    filtered_graph, pos = na.calculate_positions_2d(
        attribute=attribute,
        scale=scale,
        threshold=threshold,
        mode=mode,
        type_relationship=type_relationship,
    )

    # Print some analysis results
    print("\nNetwork Analysis Results:")
    print(f"Number of nodes: {filtered_graph.number_of_nodes()}")
    print(f"Number of edges: {filtered_graph.number_of_edges()}")

    # Calculate some network metrics
    node_coefficients, avg_clustering = na.compute_clustering_coefficients()
    print(f"\nAverage clustering coefficient: {avg_clustering:.3f}")

    connectivity = na.analyze_connectivity()
    print(f"\nIs network connected? {connectivity['is_connected']}")
    print(f"Number of components: {len(connectivity['components'])}")

    # Analyze a specific node (first one in the graph)
    if filtered_graph.nodes:
        first_node = list(filtered_graph.nodes)[0]
        node_analysis = na.analyze_node(first_node)
        print(f"\nAnalysis of node {first_node}:")
        print(f"Degree: {node_analysis['degree']}")
        print(f"Number of neighbors: {len(node_analysis['neighbors'])}")

    # detect communities
    communities = na.detect_communities()
    print(f"Number of communities: {len(communities)}")
    for i, community in enumerate(communities):
        print(f"Community {i+1}: {community}")
