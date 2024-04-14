import pandas as pd
from tqdm import tqdm
import networkx as nx
import plotly.express as px
import plotly.graph_objects as go
from typing import List, Optional
from itertools import combinations
from pydantic import BaseModel, Field
from Bio.Align import Alignment as BioAlignment
from joblib import Parallel, delayed, cpu_count
from typing import List, Optional, Union, Tuple, TYPE_CHECKING
from concurrent.futures import ThreadPoolExecutor, as_completed


from pyeed.core.sequence import Sequence
from pyeed.core.abstractsequence import AbstractSequence
from pyeed.core.pairwisealignment import PairwiseAlignment
from pyeed.core.proteininfo import ProteinInfo
from pyeed.aligners.pairwise import PairwiseAligner


class SequenceNetwork(BaseModel):
    """
    A class representing a sequence network.

    The SequenceNetwork class is used to create and visualize a network of sequences. It takes a list of AbstractSequence objects. The network can be visualized in 2D or 3D, with nodes representing sequences and edges representing alignments between sequences.

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

    targets: Optional[List[str]] = Field(
        default=[],
        description="List of selected sequences",
    )

    def __init__(self, sequences: List[AbstractSequence], weight: str = "identity", color: str = "name", threshold: float = None, label: str = "name", dimensions: int = 3):
        super().__init__()
        self.sequences = sequences
        self.threshold = threshold


    def add_target(self, target: AbstractSequence):
        if target.source_id not in self.targets:
            self.targets.append(target.source_id)

    def _create_pairwise_alignments(self, input_sequences, aligner: "PairwiseAligner", **kwargs):
        """
        Creates pairwise alignments between sequences.

        This method creates pairwise alignments between sequences in the network.
        The pairwise alignments are stored in the 'pairwise_alignments' attribute of the SequenceNetwork object.
        This is done for the later visualization of the network graph with cytoscope.

        Args:
            aligner (PairwiseAligner): Python-based aligner to be called.

        Raises:
            ValueError: If the number of sequences is less than 2.

        Returns:
            Nothing the data is stored internally in fields of the class.
        """

        # Pairwise alignment
        if len(input_sequences) == 2:

            pairwise_aligner = aligner(
                sequences=[
                    input_sequences[0].sequence,
                    input_sequences[1].sequence,
                ],
                **kwargs,
            )
            alignment_result = pairwise_aligner.align()

            return self._map_pairwise_alignment_results(
                alignment_result,
                pair=(
                    input_sequences[0],
                    input_sequences[1],
                ),
                mode=pairwise_aligner.mode,
            )

        # Multi pairwise alignment
        elif len(input_sequences) > 2:
            pairs = list(combinations(input_sequences, 2))

            aligners = [
                aligner(sequences=[s.sequence for s in pair], **kwargs)
                for pair in pairs
            ]

            alignments = Parallel(n_jobs=cpu_count(), prefer="processes")(
                delayed(a.align)()
                for a in tqdm(aligners, desc="⛓️ Running pairwise alignments")
            )

            return alignments, pairs, aligners[0].mode
                

        else:
            raise ValueError(
                f"Alignment Error. Recieved {len(input_sequences)} sequences. Expected 2."
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
        alignments, pairs, mode = self._create_pairwise_alignments(self.sequences, PairwiseAligner, mode="global")

        # Add nodes and assign node attributes
        # TODO die Aminosäuren sequenz soll jetzt auch noch da rein
        if all([isinstance(sequence, ProteinInfo) for sequence in self.sequences]):
            
            node_data = []
            
            def process_sequence(sequence):
                return (sequence.source_id, {
                    'name': sequence.name,
                    'family_name': sequence.family_name,
                    'domain': sequence.organism.domain,
                    'kingdom': sequence.organism.kingdom,
                    'phylum': sequence.organism.phylum,
                    'tax_class': sequence.organism.tax_class,
                    'order': sequence.organism.order,
                    'family': sequence.organism.family,
                    'genus': sequence.organism.genus,
                    'species': sequence.organism.species,
                    'ec_number': sequence.ec_number,
                    'mol_weight': sequence.mol_weight,
                    'taxonomy_id': sequence.organism.taxonomy_id,
                })

            # Setting up a ThreadPoolExecutor to manage a pool of threads
            with ThreadPoolExecutor(max_workers=10) as executor:
                # Submitting tasks to the executor
                futures = [executor.submit(process_sequence, seq) for seq in self.sequences]

                # Collecting results as they complete
                for future in as_completed(futures):
                    node_data.append(future.result())



            # for sequence in self.sequences:

            #     node_data.append((sequence.source_id, {
            #         'name': sequence.name,
            #         'family_name': sequence.family_name,
            #         'domain': sequence.organism.domain,
            #         'kingdom': sequence.organism.kingdom,
            #         'phylum': sequence.organism.phylum,
            #         'tax_class': sequence.organism.tax_class,
            #         'order': sequence.organism.order,
            #         'family': sequence.organism.family,
            #         'genus': sequence.organism.genus,
            #         'species': sequence.organism.species,
            #         'ec_number': sequence.ec_number,
            #         'mol_weight': sequence.mol_weight,
            #         'taxonomy_id': sequence.organism.taxonomy_id,
            #     }))                

                # graph.add_node(
                #     sequence.source_id,
                #     name=sequence.name,
                #     familiy_name=sequence.family_name,
                #     domain=sequence.organism.domain,
                #     kingdome=sequence.organism.kingdom,
                #     phylum=sequence.organism.phylum,
                #     tax_class=sequence.organism.tax_class,
                #     order=sequence.organism.order,
                #     family=sequence.organism.family,
                #     genus=sequence.organism.genus,
                #     species=sequence.organism.species,
                #     ec_number=sequence.ec_number,
                #     mol_weight=sequence.mol_weight,
                #     taxonomy_id=sequence.organism.taxonomy_id,
                # )

                graph.add_nodes_from(node_data)



        else:
            for sequence in self.sequences:
                graph.add_node(
                    sequence.source_id,
                    name=sequence.name,
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

        print('this is the threshold', self.threshold)


        if self.threshold != None:

            # create a datafarme for the egdes
            edge_data = []

            def process_pair(alignment, pair):
                shorter_seq = min(pair, key=lambda x: len(x.sequence))
                
                identities = alignment.counts().identities
                identity = identities / len(shorter_seq.sequence)
                
                if identity >= self.threshold:
                    return (
                        pair[0].source_id, pair[1].source_id,
                        {
                            'identity': identity,
                            'gaps': 1 / (alignment.counts().gaps + 1),
                            'mismatches': 1 / (alignment.counts().mismatches + 1),
                            'score': alignment.score
                        }
                    )
                return None

            with ThreadPoolExecutor() as executor:
                # Submit tasks to the executor.
                future_to_pair = {executor.submit(process_pair, alignment, pair): (alignment, pair) for alignment, pair in zip(alignments, pairs)}
                
                # Collecting results as they complete.
                for future in as_completed(future_to_pair):
                    result = future.result()
                    if result is not None:
                        edge_data.append(result)



            # for alignment, pair in zip(alignments, pairs):
            #     shorter_seq = min(list(pair), key=lambda x: len(x.sequence))

            #     identities = alignment.counts().identities
            #     identity = identities / len(shorter_seq.sequence)

            #     if identity >= self.threshold:
            #         edge_data.append(((list(pair)[0].source_id, list(pair)[1].source_id,{ 
            #             'identity': identity, 
            #             'gaps': 1 / (alignment.counts().gaps + 1), 
            #             'mismatches': 1 / (alignment.counts().mismatches + 1), 
            #             'score': alignment.score} )))
            
            graph.add_edges_from(edge_data)
                    

        else:
            for alignment, pair in zip(alignments, pairs):

                identities = alignment.counts().identities
                identity = identities / len(shorter_seq.sequence)

                graph.add_edge(
                        list(pair)[0].source_id,
                        list(pair)[1].source_id,
                        identity=identity,
                        gaps = 1 / (alignment.counts().gaps + 1),
                        mismatches = 1 / (alignment.counts().mismatches + 1),
                        score = alignment.score,
                )

        # Calculate betweenness centrality
        betweenness = nx.betweenness_centrality(graph)
        for node, betw_value in betweenness.items():
            graph.nodes[node]["betweenness"] = betw_value

        # Calculate node positions based on dimensions
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
