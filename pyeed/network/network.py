import time
from typing import List, Optional

import networkx as nx
import py4cytoscape as p4c
import plotly.express as px
import plotly.graph_objects as go
from pydantic import BaseModel, Field, PrivateAttr

from pyeed.core.proteinrecord import ProteinRecord
from pyeed.core.sequencerecord import SequenceRecord
from pyeed.align.pairwise_aligner import PairwiseAligner


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

    sequences: Optional[List[SequenceRecord]] = Field(
        default=[],
        description="List of sequences to be compared",
    )

    mode: Optional[str] = Field(
        default="global",
        description="Alignment mode",
    )

    weight: Optional[str] = Field(
        default="identity",
        description="Attribute of Alignment to weight the edges",
    )

    dimensions: Optional[int] = Field(
        default=3,
        description="Dimension of the network graph",
    )

    targets: Optional[List[str]] = Field(
        default=[],
        description="List of selected sequences",
    )

    network: Optional[nx.Graph] = Field(
        default=None,
        description="Network graph with networkx",
    )

    _aligner: Optional[PairwiseAligner] = PrivateAttr(
        default=None,
    )

    def __init__(
        self,
        sequences: List[SequenceRecord],
        weight: str = "identity",
        dimensions: int = 3,
        mode = "global"
    ):
        super().__init__()
        self.weight = weight
        self.dimensions = dimensions
        self.mode = mode
        self.sequences = sequences
        self.network = nx.Graph()
        self._aligner = PairwiseAligner(mode=self.mode)

        self._create_graph()


    def add_target(self, target: SequenceRecord):
        # TODO find out what to do with targets
        if target.source_id not in self.targets:
            self.targets.append(target.source_id)

    def _process_sequence(self, sequence: ProteinRecord):
        return (
            sequence.id,
            sequence.sequence,
            {
                "name": sequence.name,
                "domain": sequence.organism.domain,
                "kingdom": sequence.organism.kingdom,
                "phylum": sequence.organism.phylum,
                "tax_class": sequence.organism.tax_class,
                "order": sequence.organism.order,
                "family": sequence.organism.family,
                "genus": sequence.organism.genus,
                "species": sequence.organism.species,
                "ec_number": sequence.ec_number,
                "mol_weight": sequence.mol_weight,
                "taxonomy_id": sequence.organism.taxonomy_id,
            },
        )

    def _create_graph(self):
        # first we add the nodes to the network
        # in the same loop we read out the sequences and the key in order to be able to perform the alignment next
        alignment_data = {}

        if all([isinstance(sequence, ProteinRecord) for sequence in self.sequences]):
            node_data = []

            print(time.time(), '----- Processing in for loop')

            for sequence in self.sequences:
                id_seq, seq, data = self._process_sequence(sequence)
                node_data.append((id_seq, data))
                alignment_data[id_seq] = seq

            print(time.time(), "----- Adding in networx")

            self.network.add_nodes_from(node_data)

            print(time.time(), "----- nodes added ")

        else:
            for sequence in self.sequences:

                id_seq = sequence.id
                seq = sequence.sequence
                alignment_data[id_seq] = seq

                self.network.add_node(
                    id_seq,
                    name=sequence.name,
                    sequence=seq,
                    domain=sequence.organism.domain,
                    kingdome=sequence.organism.kingdom,
                    phylum=sequence.organism.phylum,
                    tax_class=sequence.organism.tax_class,
                    order=sequence.organism.order,
                    family=sequence.organism.family,
                    genus=sequence.organism.genus,
                    species=sequence.organism.species,
                    taxonomy_id=sequence.organism.taxonomy_id,
                )

        # create the alignments
        alignments_results = self._aligner.align_multipairwise(alignment_data)
        # now pydantic cant seem to handle the serialization of the alignment results in the field start and end
        # so we need to convert them to lists
        for alignment_result in alignments_results:
            alignment_result['start'] = alignment_result['start'].tolist()
            alignment_result['end'] = alignment_result['end'].tolist()
        # create a list for the egdes
        edge_data = []

        print(time.time(), "Processing Pairs in for loop")
        
        for alignment_result in alignments_results:
            edge = (alignment_result['seq1_id'], alignment_result['seq2_id'], {key: value for key, value in alignment_result.items() if key not in ['seq1_id', 'seq2_id']})
            if edge:
                edge_data.append(edge)
    
        print(time.time(), "Adding edeges from list in networx")

        self.network.add_edges_from(edge_data)

        print(time.time(), "edges added in network")
        # Calculate betweenness centrality
        betweenness = nx.betweenness_centrality(self.network)
        for node, betw_value in betweenness.items():
            self.network.nodes[node]["betweenness"] = betw_value

        # Calculate degree of nodes without filtering
        nx.set_node_attributes(self.network, dict(nx.degree(self.network)), "degree_all")

        # Calculate node positions based on dimensions
        if self.dimensions == 2:
            return self._2d_position_nodes_and_edges(self.network)
        elif self.dimensions == 3:
            return self._3d_position_nodes_and_edges(self.network)
        else:
            if self.dimensions > 3:
                raise ValueError(
                    f"Bruuuhh chill, u visiting from {self.dimensions}D cyberspace? Dimensions must be 2 or 3"
                )
            
    def create_cytoscape_graph(
        self, collection: str, title: str, threshold: float = 0.8
    ):
        # assert that the cytoscape API is running and cytoscape is running in the background
        assert p4c.cytoscape_ping(), "Cytoscape is not running in the background"
        assert p4c.cytoscape_version_info(), "Cytoscape API is not running"
        # filter the the edges by the threshold


        p4c.create_network_from_networkx(
            self.network, collection=collection, title=title + "_" + str(threshold)
        )

    def set_layout(
        self,
        layout_name: str = "force-directed",
        properties_dict: dict = {
            "defaultSpringCoefficient": 4e-5,
            "defaultSpringLength": 100,
            "defaultNodeMass": 3,
            "numIterations": 50,
        },
    ):
        p4c.layout_network(layout_name)
        # ['numIterations', 'defaultSpringCoefficient', 'defaultSpringLength', 'defaultNodeMass', 'isDeterministic', 'singlePartition']
        p4c.set_layout_properties(
            layout_name=layout_name, properties_dict=properties_dict
        )

        p4c.scale_layout(axis="Both Axis", scale_factor=1.0)
        p4c.layout_network(layout_name)

    def filter_cytoscape_edges_by_parameter(self, name: str, parameter: str, min_val: float, max_val: float):
        # this is a filter for the network in cytoscape
        # here the nodes and edges not relevant are filtered out
        p4c.create_column_filter(
            name,
            parameter,
            [min_val, max_val],
            "BETWEEN",
            type="edges",
            apply=True,
            hide=True,
        )

    def calculate_degree(self):
        # Calculate degree of nodes with filtering
        g = p4c.create_networkx_from_network()
        nx.set_node_attributes(g, dict(nx.degree(self.network)), "degree")
        p4c.create_network_from_networkx(g, collection="tet", title="hfjakd")

    def visualize(self):
        """
        Visualizes the network graph.

        This method visualizes the network graph created by the SequenceNetwork class.
        It checks the value of the 'dimensions' attribute and calls either the 'visualize_2d_network' or 'visualize_3d_network' method accordingly.
        If the 'dimensions' attribute is not 2 or 3, it raises a ValueError.

        Returns:
            None

        Raises:
            ValueError: If the 'dimensions' attribute is greater than 3.
        """

        if self.dimensions == 2:
            self.visualize_2d_network()
        elif self.dimensions == 3:
            self.visusalize_3d_network()
        else:
            if self.dimensions > 3:
                raise ValueError(
                    f"Dimensions must be 2 or 3 for visualization, not {self.dimensions}"
                )

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

    def visusalize_3d_network(self):
        """Visualizes a 3D network graph."""

        traces = []

        betweenness = 1
        bridge_nodes = [
            node_id
            for node_id, betw in self.graph.nodes(data="betweenness")
            if betw > betweenness
        ]
        for edge in self.graph.edges:
            if edge[0] in bridge_nodes or edge[1] in bridge_nodes:
                traces.append(
                    go.Scatter3d(
                        x=self.graph.edges[edge]["x_pos"],
                        y=self.graph.edges[edge]["y_pos"],
                        z=self.graph.edges[edge]["z_pos"],
                        mode="lines",
                        hoverinfo="none",
                        line=dict(
                            width=1,
                            color="rgba(128, 128, 128, 0.8)",
                        ),
                    )
                )

        color_conditions = set([node[self.color] for node in self.graph.nodes.values()])
        color_dict = dict(
            zip(color_conditions, self._sample_colorscale(len(color_conditions)))
        )

        # Add nodes
        for key, node in self.graph.nodes.items():
            size = 6
            color = color_dict[node[self.color]]
            symbol = "circle"
            if key in self.targets:
                size = 10
                color = "red"
                symbol = "cross"

            info = [tuple(n for n in node.values())]

            traces.append(
                go.Scatter3d(
                    x=[node["x_pos"]],
                    y=[node["y_pos"]],
                    z=[node["z_pos"]],
                    mode="markers",
                    marker=dict(
                        size=size,
                        color=color,
                        symbol=symbol,
                    ),
                    text=node["name"],
                    hovertemplate="<b>%{customdata[0]}</b><br>Family Name: %{customdata[1]}<br>Domain: %{customdata[2]}<br>Kingdom: %{customdata[3]}</b><br>Phylum: %{customdata[4]}<br>Class: %{customdata[5]}<br>Order: %{customdata[6]}<br>Family: %{customdata[7]}<br>Genus: %{customdata[8]}<br>Species: %{customdata[9]}<br>EC Number: %{customdata[10]}<br>Mol Weight: %{customdata[11]}<br>Taxonomy ID: %{customdata[12]}<extra></extra>",
                    customdata=list((info)),
                )
            )

        # Plot figure
        fig = go.Figure(
            data=traces,
            layout=go.Layout(
                plot_bgcolor="white",
                showlegend=False,
                hovermode="closest",
                margin=dict(b=0, l=0, r=0, t=0),
            ),
        )
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)
        fig.update_scenes(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False)

        fig.show()

    def visualize_2d_network(self):
        """Visualizes a 2D network graph."""
        traces = []
        betweenness = 1
        bridge_nodes = [
            node_id
            for node_id, betw in self.graph.nodes(data="betweenness")
            if betw > betweenness
        ]
        for edge in self.graph.edges:
            if edge[0] in bridge_nodes or edge[1] in bridge_nodes:
                traces.append(
                    go.Scatter(
                        x=self.graph.edges[edge]["x_pos"],
                        y=self.graph.edges[edge]["y_pos"],
                        mode="lines",
                        line=dict(
                            width=1,
                            color="rgba(128, 128, 128, 0.1)",
                        ),
                    )
                )

        # Add nodes
        color_conditions = set([node[self.color] for node in self.graph.nodes.values()])
        color_dict = dict(
            zip(color_conditions, self._sample_colorscale(len(color_conditions)))
        )

        for key, node in self.graph.nodes.items():
            size = 6
            color = color_dict[node[self.color]]
            symbol = "circle"
            if key in self.targets:
                size = 10
                color = "red"
                symbol = "cross"

            info = [tuple(n for n in node.values())]
            traces.append(
                go.Scatter(
                    x=[node["x_pos"]],
                    y=[node["y_pos"]],
                    mode="markers",
                    marker=dict(
                        size=size,
                        color=color,
                        symbol=symbol,
                    ),
                    text=node["name"],
                    hovertemplate="<b>%{customdata[0]}</b><br>Family Name: %{customdata[1]}<br>Domain: %{customdata[2]}<br>Kingdom: %{customdata[3]}</b><br>Phylum: %{customdata[4]}<br>Class: %{customdata[5]}<br>Order: %{customdata[6]}<br>Family: %{customdata[7]}<br>Genus: %{customdata[8]}<br>Species: %{customdata[9]}<br>EC Number: %{customdata[10]}<br>Mol Weight: %{customdata[11]}<br>Taxonomy ID: %{customdata[12]}<extra></extra>",
                    customdata=list((info)),
                )
            )

        # Plot figure
        fig = go.Figure(
            data=traces,
            layout=go.Layout(
                plot_bgcolor="white",
                showlegend=False,
                hovermode="closest",
                margin=dict(b=0, l=0, r=0, t=0),
            ),
        )
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)

        fig.show(scale=10)

    @staticmethod
    def _sample_colorscale(size: int) -> List[str]:
        return px.colors.sample_colorscale("viridis", [i / size for i in range(size)])
