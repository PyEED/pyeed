from typing import List, Optional
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from pydantic import BaseModel, Field

from pyeed.core.abstractsequence import AbstractSequence
from pyeed.core.pairwisealignment import PairwiseAlignment


class SequenceNetwork(BaseModel):
    """
    A class representing a sequence network.

    The SequenceNetwork class is used to create and visualize a network of sequences. It takes a list of AbstractSequence objects and a list of PairwiseAlignment objects as input. The network can be visualized in 2D or 3D, with nodes representing sequences and edges representing alignments between sequences.

    Attributes:
        sequences (Optional[List[AbstractSequence]]): A list of AbstractSequence objects to be compared in the network. Default is an empty list.
        pairwise_alignments (Optional[List[PairwiseAlignment]]): A list of PairwiseAlignment objects representing the pairwise alignments between sequences. Default is an empty list.
        weight (Optional[str]): The attribute of the Alignment object to weight the edges in the network. Default is "identity".
        color (Optional[str]): The attribute of the ProteinInfo object to colorize the nodes in the network. Default is "name".
        threshold (Optional[float]): Sequences with a weight higher than the threshold are connected in the network. Default is None.
        label (Optional[str]): The node label in the graph. Default is "name".
        dimensions (Optional[int]): The dimension of the network graph. Default is 3.

    Methods:
        graph() -> nx.Graph: Maps properties of alignments to a network graph.
        visualize(): Visualizes the network graph.
    """

    sequences: Optional[List[AbstractSequence]] = Field(
        default=[],
        description="List of sequences to be compared",
    )

    pairwise_alignments: Optional[List[PairwiseAlignment]] = Field(
        default=[],
        description="List of pairwise alignments",
    )

    weight: Optional[str] = Field(
        default="identity",
        description="Attribute of Alignment to weight the edges",
    )

    color: Optional[str] = Field(
        default="name",
        description="Attribute of ProteinInfo to colorize nodes",
    )

    threshold: Optional[float] = Field(
        default=None,
        description="Sequences with a weight higher than the threshold are connected in the network",
    )

    label: Optional[str] = Field(
        default="name",
        description="Node label in the graph",
    )

    dimensions: Optional[int] = Field(
        default=3,
        description="Dimension of the network graph",
    )

    @property
    def graph(self) -> nx.Graph:
        """
        Maps properties of alignments to a network graph.

        Returns:
            nx.Graph: The network graph representing the sequence network.

        Raises:
            ValueError: If the dimensions of the network graph are greater than 3.

        Notes:
            - The graph is created using the NetworkX library.
            - Nodes in the graph represent sequences, and edges represent alignments between sequences.
            - Node attributes include the sequence name, organism, and taxonomy ID.
            - Edge attributes include the alignment identity, gaps, mismatches, and score.
            - The graph can be visualized in 2D or 3D using the visualize() method.
        """
        graph = nx.Graph()

        # Add nodes and assign node attributes
        for sequence in self.sequences:
            graph.add_node(
                sequence.source_id,
                name=sequence.name,
                organism=sequence.organism,
                taxonomy_id=sequence.organism.taxonomy_id,
            )

        # Add edges and assign edge attributes
        if self.threshold != None:
            for alignment in self.pairwise_alignments:
                if alignment.identity >= self.threshold:
                    graph.add_edge(
                        alignment.input_sequences[0].source_id,
                        alignment.input_sequences[1].source_id,
                        identity=alignment.identity,
                        gaps=1 / (alignment.gaps + 1),
                        mismatches=1 / (alignment.mismatches + 1),
                        score=alignment.score,
                    )
        else:
            for alignment in self.pairwise_alignments:
                graph.add_edge(
                    alignment.input_sequences[0].source_id,
                    alignment.input_sequences[1].source_id,
                    identity=alignment.identity,
                    gaps=1 / (alignment.gaps + 1),
                    mismatches=1 / (alignment.mismatches + 1),
                    score=alignment.score,
                )

        if self.dimensions == 2:
            return self._2d_position_nodes_and_edges(graph)
        elif self.dimensions == 3:
            return self._3d_position_nodes_and_edges(graph)
        else:
            if self.dimensions > 3:
                raise ValueError(
                    f"Bruuuhh chill, u visiting from {self.dimensions}D cyberspace? Dimensions must be 2 or 3"
                )

        return graph

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

        Notes:
            - The visualization is done using the Plotly library.
            - If the 'dimensions' attribute is 2, the network graph is visualized in 2D.
            - If the 'dimensions' attribute is 3, the network graph is visualized in 3D.
            - The visualization includes nodes representing sequences and edges representing alignments between sequences.
            - The color of the nodes is determined by the 'color' attribute of the SequenceNetwork class.
            - The hover information for each node includes the sequence name, organism, kingdom, phylum, and source ID.
            - The visualization is displayed using the Plotly library.
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
        # Add edges
        # for edge in self.graph.edges:
        #     traces.append(
        #         go.Scatter3d(
        #             x=self.graph.edges[edge]["x_pos"],
        #             y=self.graph.edges[edge]["y_pos"],
        #             z=self.graph.edges[edge]["z_pos"],
        #             mode="lines",
        #             line=dict(
        #                 width=1,
        #                 color="rgba(128, 128, 128, 0.1)",
        #             ),
        #             hoverinfo="text",
        #             text=self.graph.edges[edge]["gaps"],
        #         )
        #     )

        color_conditions = set([node[self.color] for node in self.graph.nodes.values()])
        color_dict = dict(
            zip(color_conditions, self._sample_colorscale(len(color_conditions)))
        )

        # Add nodes
        for key, node in self.graph.nodes.items():
            info = [
                (
                    node["name"],
                    node["organism"].name,
                    node["organism"].kingdom,
                    node["organism"].phylum,
                    key,
                )
            ]

            traces.append(
                go.Scatter3d(
                    x=[node["x_pos"]],
                    y=[node["y_pos"]],
                    z=[node["z_pos"]],
                    mode="markers",
                    marker=dict(
                        size=7,
                        color=color_dict[node[self.color]],
                    ),
                    text=node["name"],
                    hovertemplate="<b>%{customdata[0]}</b><br>Organism: %{customdata[1]}<br>Kingdom: %{customdata[2]}<br>Phylum: %{customdata[3]}</b><br>Source ID: %{customdata[4]}<extra></extra>",
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

        # Add nodes
        color_conditions = set([node[self.color] for node in self.graph.nodes.values()])
        color_dict = dict(
            zip(color_conditions, self._sample_colorscale(len(color_conditions)))
        )

        for key, node in self.graph.nodes.items():
            info = [
                (
                    node["name"],
                    node["organism"].name,
                    node["organism"].kingdom,
                    node["organism"].phylum,
                    key,
                )
            ]

            traces.append(
                go.Scatter(
                    x=[node["x_pos"]],
                    y=[node["y_pos"]],
                    mode="markers",
                    marker=dict(
                        size=10,
                        color=color_dict[node[self.color]],
                    ),
                    text=node["name"],
                    hovertemplate="<b>%{customdata[0]}</b><br>Organism: %{customdata[1]}<br>Kingdom: %{customdata[2]}<br>Phylum: %{customdata[3]}</b><br>Source ID: %{customdata[4]}<extra></extra>",
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

        fig.show()

    @staticmethod
    def _sample_colorscale(size: int) -> List[str]:
        return px.colors.sample_colorscale("viridis", [i / size for i in range(size)])
