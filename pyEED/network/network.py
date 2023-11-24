from typing import List
import networkx as nx
import plotly.graph_objects as go

from pyEED.core.proteininfo import ProteinInfo
from pyEED.core.pairwisealignment import PairwiseAlignment


def pairwise_network(
    alignments: List[PairwiseAlignment],
    weight: str = "identity",
    cutoff: float = None,
    label: str = "accession_id",
    color: str = "name",
) -> None:
    """Takes a list of `PairwiseAlignment`, constructs a network graph,
    whereas the edges of the graph are weighted with an attribute of the
    `PairwiseAlignment` object. Visualizes the network graph.

    Args:
        alignments (List[PairwiseAlignment]): List of pairwise alignments.
        weight (str, optional): Attribute of `PairwiseAlignment` to weight the edges.
        Defaults to "identity".
        label (str, optional): Node labe in the graph. Defaults to "accession_id".
        color (str, optional): Attribute of `ProteinInfo` co colorize nodes. Defaults to "name".

    Raises:
        ValueError: If 'weights' is not one of "identity", "score", "gaps", or "mismatches".
        ValueError: If 'label' is not one of "name", "organism", "ec_number",
        "mol_weight", or "accession_id".
    """

    # Validate weights
    weights = ["identity", "score", "gaps", "mismatches"]
    if weight not in weights:
        raise ValueError(
            f"'weight' to parameterize network must be an alignment property ({weights})"
        )

    # Validate labels
    labels = ["name", "organism", "ec_number", "mol_weight", "accession_id"]
    if label not in labels:
        raise ValueError(
            f"'label' to parameterize network must be an alignment property ({labels})"
        )

    graph = construct_network(alignments, cutoff=cutoff)

    graph = position_nodes_and_edges(graph, weight=weight)

    visualize_network(graph, label=label, color=color)


def construct_network(alignments: List[PairwiseAlignment], cutoff: float) -> nx.Graph:
    """Maps properties of alignments to a network graph."""
    graph = nx.Graph()

    sequences = _get_unique_sequences(alignments)

    # Add nodes and assign node attributes
    for sequence in sequences:
        graph.add_node(
            node_for_adding=sequence.source_id,
            name=sequence.name,
            accession_id=sequence.source_id,
            organism=sequence.organism.name,
            ec_number=sequence.ec_number,
            mol_weight=sequence.mol_weight,
        )

    # Add edges and assign edge attributes
    if cutoff != None: 

        for alignment in alignments:

            if alignment.identity >= cutoff:
                
                graph.add_edge(
                    alignment.reference_seq.source_id,
                    alignment.query_seq.source_id,
                    identity=alignment.identity,
                    gaps=1 / (alignment.gaps + 1),
                    mismatches=1 / (alignment.mismatches + 1),
                    score=alignment.score,
                )
    else:
        for alignment in alignments:

            graph.add_edge(
                alignment.reference_seq.source_id,
                alignment.query_seq.source_id,
                identity=alignment.identity,
                gaps=1 / (alignment.gaps + 1),
                mismatches=1 / (alignment.mismatches + 1),
                score=alignment.score,
            )

    return graph


def visualize_network(graph: nx.Graph, label: str, color: str):
    """Visualizes a network graph."""
    # TODO: Refactor coloration of nodes
    # TODO: Add ability to colorize edges

    colors = ["blue", "red", "green", "orange", "purple", "pink", "brown", "yellow"]
    traces = []
    # Add edges
    for edge in graph.edges:
        traces.append(
            go.Scatter(
                x=graph.edges[edge]["x_pos"],
                y=graph.edges[edge]["y_pos"],
                mode="lines",
                line=dict(
                    width=1,
                    # colorscale="YlGnBu",
                    # colorbar=dict(title='Heatmap Colorscale')
                    color="rgba(128, 128, 128, 0.1)",
                ),
                hoverinfo="text",
                text=graph.edges[edge]["gaps"],
            )
        )

    # Add nodes
    color_dict = {}
    for _, data in graph.nodes(data=True):
        if data[color] not in color_dict:
            color_dict[data[color]] = colors[len(color_dict)]

        traces.append(
            go.Scatter(
                x=[data["x_pos"]],
                y=[data["y_pos"]],
                mode="markers",
                hoverinfo="text",
                marker=dict(
                    # showscale=True,
                    # colorscale="Viridis",
                    size=10,
                    # colorbar=dict(thickness=15, title=""),
                    color=color_dict[data[color]],
                ),
                text=data[label],
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
    print(color_dict)

    fig.show()


def position_nodes_and_edges(graph: nx.Graph, weight: str) -> nx.Graph:
    """Calculates node positions based on weight metric and
    adds position information of nodes and edges to the respective
    entry in the graph."""

    positions = nx.spring_layout(graph, weight=weight)

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


def _get_unique_sequences(alignments: List[PairwiseAlignment]) -> List[ProteinInfo]:
    """Gets unique sequences from a list of alignments."""

    sequence_infos = []
    added_ids = set()
    for alignment in alignments:
        if alignment.reference_seq.source_id not in added_ids:
            sequence_infos.append(alignment.reference_seq)
            added_ids.add(alignment.reference_seq.source_id)

        if alignment.query_seq.source_id not in added_ids:
            sequence_infos.append(alignment.query_seq)
            added_ids.add(alignment.query_seq.source_id)

    return sequence_infos
