import networkx as nx
from pyeed.dbconnect import DatabaseConnector
from loguru import logger

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

    def create_graph(self, node_types=None, relationships=None, ids=None):
        """
        Creates a graph using data from the Neo4j database with specified filters.

        Args:
            node_types (list[str], optional): List of node types to include (e.g., ['Protein', 'DNA']).
            relationships (list[str], optional): List of relationship types to include (e.g., ['HAS_STANDARD_NUMBERING']).
            ids (list[str], optional): List of accession IDs to include specific nodes.

        Returns:
            networkx.Graph: The created graph.
        """

        logger.info(f"Creating graph with node types: {node_types} and relationships: {relationships} and ids: {ids}")

        # Query to fetch nodes with filters
        node_filter = ""
        if node_types:
            node_filter += f"WHERE labels(n)[0] IN $node_types "
        if ids:
            if "WHERE" in node_filter:
                node_filter += f"AND n.accession_id IN $accession_ids "
            else:
                node_filter += f"WHERE n.accession_id IN $accession_ids "

        query_nodes = f"""
        MATCH (n)
        {node_filter}
        RETURN ID(n) as id, labels(n) as labels, properties(n) as properties
        """

        # Query to fetch relationships with filters
        relationship_filter = ""
        if relationships:
            relationship_filter += f"WHERE type(r) IN $relationships "

        query_relationships = f"""
        MATCH (n)-[r]->(m)
        {relationship_filter}
        RETURN ID(n) as source, ID(m) as target, type(r) as type, properties(r) as properties
        """

        # Fetch nodes and relationships
        logger.debug(f"Executing query: {query_nodes}")
        nodes = self.db.execute_read(query_nodes, {"node_types": node_types, "accession_ids": ids})
        logger.debug(f"Executing query: {query_relationships}")
        relationships = self.db.execute_read(query_relationships, {"relationships": relationships})
        logger.debug(f"Number of nodes: {len(nodes)}")
        logger.debug(f"Number of relationships: {len(relationships)}")

        # Add nodes
        for node in nodes:
            self.graph.add_node(
                node["id"],
                labels=node["labels"],
                properties=node["properties"],
            )

        # Add relationships
        for rel in relationships:
            if rel["source"] in self.graph and rel["target"] in self.graph:
                self.graph.add_edge(
                    rel["source"],
                    rel["target"],
                    type=rel["type"],
                    properties=rel["properties"],
                )

        return self.graph

    def compute_degree_centrality(self):
        """
        Computes the degree centrality of the graph.

        Returns:
            dict: A dictionary mapping each node to its degree centrality.
        """
        return nx.degree_centrality(self.graph)

    def shortest_path(self, source, target):
        """
        Finds the shortest path between two nodes.

        Args:
            source (node): The source node.
            target (node): The target node.

        Returns:
            list: A list of nodes representing the shortest path.
        """
        try:
            return nx.shortest_path(self.graph, source=source, target=target)
        except nx.NetworkXNoPath:
            return None

    def detect_communities(self):
        """
        Detects communities in the graph using the greedy modularity algorithm.

        Returns:
            list: A list of sets, each containing the nodes in a community.
        """
        from networkx.algorithms.community import greedy_modularity_communities

        return list(greedy_modularity_communities(self.graph))

    def find_isolated_nodes(self, graph) -> list:
        """
        Finds isolated nodes in the graph.

        Args:
            graph (networkx.Graph): The graph to analyze.

        Returns:
            list: List of nodes that have no connections (degree = 0).
        """
        return [node for node in graph if graph.degree(node) == 0]
    
    def find_self_referential_nodes(self, relationship_type: str, graph: nx.Graph) -> list:
        """
        Finds all edges in the graph that are self-referential (loops) of a specified relationship type.

        Args:
            relationship_type (str): The relationship type to consider.
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
        scale: float = 1.0,
        threshold: float | None = None,
        path: str | None = None,
        mode: str = "HIDE_UNDER_THRESHOLD",
        type_relationship: str | None = None
    ) -> tuple[nx.Graph, dict[str, tuple[float, float]]]:
        """
        Calculates 2D positions for graph visualization with optional edge filtering.

        Args:
            attribute (str, optional): Edge attribute to use for weighting.
            scale (float): Scale factor for node distances in the visualization. Defaults to 1.0.
            threshold (float, optional): Threshold value for edge filtering.
            path (str, optional): Path to save the visualization (not used in current implementation).
            mode (str): Filtering mode - either "HIDE_UNDER_THRESHOLD" or "HIDE_OVER_THRESHOLD". Defaults to "HIDE_UNDER_THRESHOLD".
            type_relationship (str, optional): Specific relationship type to filter on.

        Returns:
            tuple[networkx.Graph, dict]: A tuple containing:
                - The filtered graph
                - Dictionary mapping node IDs to their 2D positions (x, y)
        """

        # Filter edges based on the threshold
        edges = self.graph.edges(data=True)
        edges_filtered = []

        if threshold is not None and attribute is not None:
            for edge in edges:
                if mode == "HIDE_UNDER_THRESHOLD":
                    if type_relationship is not None:
                        if edge[2]['type'] == type_relationship:
                            if edge[2]['properties'][attribute] >= threshold:
                                edges_filtered.append(edge)
                    else:
                        if edge[2]['properties'][attribute] >= threshold:
                            edges_filtered.append(edge)
                elif mode == "HIDE_OVER_THRESHOLD":
                    if type_relationship is not None:
                        if edge[2]['type'] == type_relationship:
                            if edge[2]['properties'][attribute] <= threshold:
                                edges_filtered.append(edge)
                    else:
                        if edge[2]['properties'][attribute] <= threshold:
                            edges_filtered.append(edge)

        logger.debug(f"Number of edges filtered: {len(edges) - len(edges_filtered)}")

        # Extract filtered graph
        filtered_graph = nx.Graph()
        filtered_graph.add_nodes_from(self.graph.nodes(data=True))
        filtered_graph.add_edges_from(edges_filtered)

        # Find self-referential nodes
        self_referential_edges = self.find_self_referential_nodes(type_relationship, filtered_graph)
        logger.info(f"Number of self-referential edges: {len(self_referential_edges)}")
        filtered_graph.remove_edges_from(self_referential_edges)

        # Find isolated nodes
        isolated_nodes = self.find_isolated_nodes(filtered_graph)
        filtered_graph.remove_nodes_from(isolated_nodes)

        # Use spring layout for force-directed graph
        pos = nx.spring_layout(filtered_graph, weight=attribute, scale=scale)

        return filtered_graph, pos
    
    def compute_clustering_coefficients(self):
        """
        Computes the clustering coefficient for all nodes and the overall graph.

        Returns:
            tuple: A tuple containing the clustering coefficients for nodes and the overall clustering coefficient.
        """
        return nx.clustering(self.graph), nx.average_clustering(self.graph)

    def analyze_connectivity(self):
        """
        Analyzes the connectivity of the graph.

        Returns:
            dict: Information about whether the graph is connected and its connected components.
        """
        is_connected = nx.is_connected(self.graph)
        components = list(nx.connected_components(self.graph)) if not is_connected else []
        return {"is_connected": is_connected, "components": components}

    def analyze_node(self, node_id):
        """
        Analyzes the properties and relationships of a specific node.

        Args:
            node_id (int): The ID of the node to analyze.

        Returns:
            dict: Information about the node.
        """
        if node_id not in self.graph:
            raise ValueError(f"Node {node_id} does not exist in the graph.")

        node_data = self.graph.nodes[node_id]
        neighbors = list(self.graph.neighbors(node_id))

        return {
            "id": node_id,
            "labels": node_data.get("labels"),
            "properties": node_data.get("properties"),
            "degree": self.graph.degree[node_id],
            "neighbors": neighbors,
        }




if __name__ == "__main__":

    uri = "bolt://localhost:7687"
    username = "neo4j"
