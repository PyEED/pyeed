import copy
from typing import List, Optional

import matplotlib.pyplot as plt
import networkx as nx
import py4cytoscape as p4c
from py4cytoscape import gen_node_size_map, scheme_c_number_continuous
from pydantic import BaseModel, Field, PrivateAttr
from requests import RequestException

from pyeed.align.pairwise import PairwiseAligner
from pyeed.core.proteinrecord import ProteinRecord
from pyeed.core.sequencerecord import SequenceRecord

plt.rcParams["figure.dpi"] = 300


class SequenceNetwork(BaseModel):
    """
    A class representing a sequence network.

    The SequenceNetwork class is used to create and visualize a network of sequences. It takes a list of SequenceRecords objects. The network can be visualized in 2D or 3D, with nodes representing sequences and edges representing alignments between sequences.

    Attributes:
        sequences (Optional[List[SequenceRecord]]): A list of AbstractSequence objects to be compared in the network. Default is an empty list.
        weight (Optional[str]): The attribute of the Alignment object to weight the edges in the network. Default is "identity".
        dimensions (Optional[int]): The dimension of the network graph. Default is 3.

    Methods:
        graph() -> nx.Graph: Maps properties of alignments to a network graph.
        visualize(): Visualizes the network graph.
    """

    class Config:
        arbitrary_types_allowed = True
        json_encoders = {
            nx.Graph: lambda g: nx.node_link_data(g),
        }

    sequences: List[SequenceRecord] = Field(
        description="List of sequences to be compared",
    )

    weight: Optional[str] = Field(
        default="identity",
        description="Attribute of Alignment to weight the edges",
    )

    targets: List[str] = Field(
        default=[],
        description="List of selected sequences",
    )

    network: Optional[nx.Graph] = Field(
        default=nx.Graph(),
        description="Network graph with networkx",
    )

    _full_network: nx.Graph = PrivateAttr(
        default=nx.Graph(),
    )

    _aligner: PairwiseAligner = PrivateAttr(
        default=PairwiseAligner(),
    )

    _cytoscape_url: Optional[str] = PrivateAttr(
        default="http://cytoscape-desktop:1234/v1",
    )

    def model_post_init(self, __context):
        self._create_graph()

    def add_to_targets(self, target: SequenceRecord):
        if target.id not in self.targets:
            self.targets.append(target.id)

    def _create_graph(self):
        """Initializes the nx.Graph object and adds nodes and edges based on the sequences."""

        sequences = {}

        # add nodes to the network
        for sequence in self.sequences:
            seq_dict = sequence.to_dict()
            seq_id = seq_dict.pop("@id")
            node_dict = seq_dict.pop("organism") if "organism" in seq_dict else {}
            node_dict["sequence"] = sequence.sequence
            node_dict["name"] = sequence.name if sequence.name else ""
            if "ec_number" in seq_dict:
                node_dict["ec_number"] = seq_dict.pop("ec_number")
            if "mol_weight" in seq_dict:
                node_dict["mol_weight"] = seq_dict.pop("mol_weight")




            sequences[seq_id] = sequence.sequence
            self._full_network.add_node(seq_id, **node_dict)

        # create the alignments
        alignments_results = self._aligner.align_multipairwise(sequences)

        # create a list for the egdes
        edge_data = []

        for alignment_result in alignments_results:
            edge = (
                alignment_result["sequences"][0]["id"],
                alignment_result["sequences"][1]["id"],
                {key: value for key, value in alignment_result.items()},
            )
            if edge:
                edge_data.append(edge)

        self._full_network.add_edges_from(edge_data)

        # Calculate node positions based on dimensions
        self._2d_position_nodes_and_edges(self._full_network)

        self.network = copy.deepcopy(self._full_network)
        self.calculate_centrality()

    def update_threshhold(self, threshold: float):
        """Removes or adds edges based on the threshold value."""

        assert 0 <= threshold <= 1, "Threshold must be between 0 and 1"

        network = copy.deepcopy(self._full_network)
        edge_pairs_below_threshold = [
            (node1, node2)
            for node1, node2, data in network.edges(data=True)
            if data["identity"] < threshold
        ]
        network.remove_edges_from(edge_pairs_below_threshold)
        self._2d_position_nodes_and_edges(network)

        self.network = network

    def calculate_centrality(self, mode: str = "betweenness"):
        """Calculates the centrality of the nodes in the network graph.


        Parameters:
            mode (str, optional): The centrality metric to calculate. Default is "betweenness".

        Returns:
            dict: A dictionary of node centrality values.
        """
        modes = ["betweenness", "closeness", "degree", "eigenvector"]

        if mode == "betweenness":
            centrality = nx.betweenness_centrality(self.network)
        elif mode == "closeness":
            centrality = nx.closeness_centrality(self.network)
        elif mode == "degree":
            centrality = nx.degree_centrality(self.network)
        elif mode == "eigenvector":
            centrality = nx.eigenvector_centrality(self.network)
        else:
            raise ValueError(
                f"Centrality mode {mode} not recognized. Choose from {modes}."
            )

        nx.set_node_attributes(self.network, centrality, f"centrality_{mode}")

    def create_cytoscape_graph(
        self,
        layout: str = "force-directed",
        threshold: float = 0.8,
        style_name: str = "default",
        column_name: str = "domain",
    ):
        try:
            p4c.cytoscape_ping(base_url=self._cytoscape_url)
        except RequestException:
            self._cytoscape_url = "http://localhost:1234/v1"
            p4c.cytoscape_ping(base_url=self._cytoscape_url)

        assert p4c.cytoscape_ping(
            base_url=self._cytoscape_url
        ), "Cytoscape is not running in the background"

        p4c.layout_network(layout, base_url=self._cytoscape_url)

        # create a degree column for the nodes based on the current chosen threshold
        self.calculate_degree(threshold=threshold)
        # filter the the edges by the threshold
        p4c.create_network_from_networkx(
            self._full_network,
            collection="SequenceNetwork",
            base_url=self._cytoscape_url,
        )

        self._hide_under_threshold(threshold)
        p4c.layout_network("grid", base_url=self._cytoscape_url)

        df_nodes = p4c.get_table_columns(table="node", base_url=self._cytoscape_url)

        data_color_names = list(set(df_nodes[column_name]))

        colors = plt.cm.tab20(range(len(data_color_names)))
        hex_colors = [
            "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))
            for r, g, b, _ in colors
        ]
        # Convert RGB to hex colors for py4cytoscape
        hex_colors = [
            "#" + "".join([f"{int(c*255):02x}" for c in color[:3]]) for color in colors
        ]

        if style_name not in p4c.get_visual_style_names(base_url=self._cytoscape_url):
            p4c.create_visual_style(style_name, base_url=self._cytoscape_url)

        p4c.set_node_color_default("#FFFFFF", style_name, base_url=self._cytoscape_url)
        p4c.set_node_color_mapping(
            column_name,
            mapping_type="discrete",
            default_color="#654321",
            style_name=style_name,
            table_column_values=data_color_names,
            colors=hex_colors,
            base_url=self._cytoscape_url,
        )

    def export_cytoscape_graph(self, file_path: str):
        """Exports the network graph to a Cytoscape-readable file.

        Args:
            file_path (str): The file path to save the network graph.

        """
        cyt_dict = nx.cytoscape_data(self.network, name="SequenceNetwork")
        with open(file_path, "w") as file:
            file.write(str(cyt_dict))
        print(f"💾 Network exported to {file_path}")

    def _hide_under_threshold(self, threshold):
        p4c.unhide_all(base_url=self._cytoscape_url)

        hide_list = []

        for u, v, d in self._full_network.edges(data=True):
            if d["identity"] < threshold:
                hide_list.append("{} (interacts with) {}".format(u, v))

        p4c.hide_edges(hide_list, base_url=self._cytoscape_url)

    def calculate_degree(self, threshold: float = 0.8):
        # Calculate degree of nodes with filtering
        degree = {}
        for u, v, d in self._full_network.edges(data=True):
            if d["identity"] > threshold:
                if u not in degree:
                    degree[u] = 1
                else:
                    degree[u] += 1
                if v not in degree:
                    degree[v] = 1
                else:
                    degree[v] += 1

        nx.set_node_attributes(
            self._full_network, degree, "degree_with_threshold_{}".format(threshold)
        )

    def set_nodes_size(
        self,
        column_name: str,
        min_size: int = 10,
        max_size: int = 100,
        style_name: str = "default",
    ):
        p4c.set_node_shape_default("ELLIPSE", style_name, base_url=self._cytoscape_url)
        # p4c.set_node_size_mapping(
        #     **gen_node_size_map(
        #         column_name,
        #         scheme_c_number_continuous(min_size, max_size),
        #         mapping_type="c",
        #         style_name=style_name,
        #         base_url=self._cytoscape_url,
        #     )
        # )
        p4c.set_node_label_mapping(
            "name", style_name=style_name, base_url=self._cytoscape_url
        )
        p4c.set_node_font_size_mapping(
            **gen_node_size_map(
                column_name,
                scheme_c_number_continuous(int(min_size / 10), int(max_size / 10)),
                style_name=style_name,
                base_url=self._cytoscape_url,
            )
        )

    def visualize(
        self,
        color: str = None,
        size: bool = False,
        edges: bool = True,
        labels: bool = False,
        save_path: str = None,
        dpi: int = 300,
    ):
        """
        Visualizes the network graph by plotting nodes and optionally edges.
        For large networks it is recommended to disable the edges for performance.

        Parameters:
            color (str, optional): The attribute to colorize the nodes. Default is None.
            size (bool, optional): Whether to size the nodes based on centrality. Default is False.
            edges (bool, optional): Whether to plot edges. Default is True.
            labels (bool, optional): Whether to add labels to the nodes. Default is False.
            save_path (str, optional): The file path to save the plot. Default is None.
            dpi (int, optional): The resolution of the saved plot. Default is 300.

        Raises:
            ValueError: If the specified color attribute is not found in the nodes.

        Returns:
            None
        """

        # plot edges
        if edges:
            [
                plt.plot(
                    edge["x_pos"],
                    edge["y_pos"],
                    color="black",
                    alpha=0.2,
                    linewidth=0.1,
                    zorder=0,
                )
                for edge in self.network.edges.values()
            ]

        # scatter positions
        node_xs = [node["x_pos"] for node in self.network.nodes.values()]
        node_ys = [node["y_pos"] for node in self.network.nodes.values()]

        # size nodes based on centrality
        if size:
            node_keys = list(self.network.nodes.values())[0].keys()
            betweenness_key = next(
                iter([key for key in node_keys if "centrality" in key])
            )
            node_sizes = [
                float((node[betweenness_key] * 300) + 1)
                for node in self.network.nodes.values()
            ]
        else:
            node_sizes = [20] * len(node_xs)

        # colorize nodes
        if color:
            try:
                color_labels = []
                for node in self.network.nodes.values():
                    color_labels.append(node[color])
            except KeyError:
                color_labels.append("n.a.")

            color_labels = list(set(color_labels))
            colors = self._sample_colorscale(len(set(color_labels)))

            color_dict = dict(
                zip(color_labels, colors)
            )

            color_list = []
            for node in self.network.nodes.values():
                try:
                    color_list.append(color_dict[node[color]])
                except KeyError:
                    color_list.append(color_dict["n.a."])

            # add legend
            markers = [
                plt.Line2D([0, 0], [0, 0], color=color, marker="o", linestyle="")
                for color in color_dict.values()
            ]
            plt.legend(
                markers,
                color_dict.keys(),
                numpoints=1,
                fontsize="xx-small",
                title=color.capitalize(),
            )
            plt.xlim(-1.1, 1.4)
        else:
            color_list = ["tab:blue"] * len(node_xs)

        if self.targets:
            for target_id in self.targets:
                target = self.network.nodes[target_id]
                plt.scatter(
                    target["x_pos"],
                    target["y_pos"],
                    c="red",
                    zorder=1,
                    marker="X",
                )

                plt.annotate(
                    target_id,
                    (target["x_pos"], target["y_pos"]),
                    fontsize=4,
                    color="black",
                )

        # plot nodes
        plt.scatter(node_xs, node_ys, c=color_list, zorder=1, s=node_sizes)

        # add labels
        if labels:
            node_names = list(self.network.nodes)
            for i, txt in enumerate(node_names):
                plt.annotate(txt, (node_xs[i], node_ys[i]), fontsize=2)

        plt.axis("off")
        if save_path:
            plt.savefig(save_path, dpi=dpi)
        plt.show()

    def _2d_position_nodes_and_edges(self, graph: nx.Graph):
        """Calculates node positions based on weight metric and
        adds position information of nodes and edges to the respective
        entry in the graph."""

        positions = nx.spring_layout(graph, weight=self.weight, dim=2, seed=42)

        # Add node position as coordinates
        for node in graph.nodes():
            graph.nodes[node]["x_pos"] = positions[node][0]
            graph.nodes[node]["y_pos"] = positions[node][1]

        # Add edge positions as coordinates
        for edge in graph.edges():
            graph.edges[edge]["x_pos"] = [
                graph.nodes[edge[0]]["x_pos"],
                graph.nodes[edge[1]]["x_pos"],
            ]
            graph.edges[edge]["y_pos"] = [
                graph.nodes[edge[0]]["y_pos"],
                graph.nodes[edge[1]]["y_pos"],
            ]

        return graph

    def _3d_position_nodes_and_edges(self, graph: nx.Graph):
        """Calculates node positions based on weight metric and
        adds position information of nodes and edges to the respective
        entry in the graph."""

        positions = nx.spring_layout(
            graph, iterations=25, weight=self.weight, dim=3, seed=42
        )

        # Add node position as coordinates
        for node in graph.nodes():
            graph.nodes[node]["x_pos"] = positions[node][0]
            graph.nodes[node]["y_pos"] = positions[node][1]
            graph.nodes[node]["z_pos"] = positions[node][2]

        # Add edge positions as coordinates
        for edge in graph.edges():
            graph.edges[edge]["x_pos"] = [
                graph.nodes[edge[0]]["x_pos"],
                graph.nodes[edge[1]]["x_pos"],
            ]
            graph.edges[edge]["y_pos"] = [
                graph.nodes[edge[0]]["y_pos"],
                graph.nodes[edge[1]]["y_pos"],
            ]
            graph.edges[edge]["z_pos"] = [
                graph.nodes[edge[0]]["z_pos"],
                graph.nodes[edge[1]]["z_pos"],
            ]

        return graph

    @staticmethod
    def _sample_colorscale(size: int):
        cmap = plt.get_cmap("viridis")
        steps = [i / size for i in range(size)]
        return [cmap(step) for step in steps]
