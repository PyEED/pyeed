# Creating Sequence Networks

A `SequenceNetwork` is created using a list of `PairwiseAlignment` objects, and a list of `AbstractSequences`. Additionally, the way the network is constructed is influenced by the `weight` and `threshold` attributes. The `weight` determines which attribute of the `PairwiseAlignment` object is used to calculate the distance between the sequences of the network. The `threshold` determines the minimum value of the `weight` attribute for an edge to be created. Furthermore, a `color` can be determined based on the attributes of an `AbstractSequence` object in which the nodes of the network will be colored.

## Visualization

=== "2D"

    ```py
    from pyeed.core import ProteinInfo, Alignment
    from pyeed.aligners import PairwiseAligner
    from pyeed.network import SequenceNetwork

    # Accessions from different methionine adenyltransferases
    mat_accessions = [
        "MBP1912539.1",
        "SEV92896.1",
        "MBO8174569.1",
        "WP_042680787.1",
        "NPA47376.1",
        "WP_167889085.1",
        "WP_048165429.1",
        "ACS90033.1",
    ]
    mats = ProteinInfo.from_ncbi(mat_accessions)

    # Create pairwise alignments between all sequences
    alignments = Alignment.from_sequences(mats, aligner=PairwiseAligner)

    # Create a network
    network = SequenceNetwork(
        sequences=mats,
        pairwise_alignments=alignments,
        weight="identity",
        threshold=0.9,
        dimensions=2,
        color="taxonomy_id",
    )

    # Visualize the network
    network.visualize()
    ```

=== "3D"

    ```py
    from pyeed.core import ProteinInfo, Alignment
    from pyeed.aligners import PairwiseAligner
    from pyeed.network import SequenceNetwork

    # Accessions from different methionine adenyltransferases
    mat_accessions = [
        "MBP1912539.1",
        "SEV92896.1",
        "MBO8174569.1",
        "WP_042680787.1",
        "NPA47376.1",
        "WP_167889085.1",
        "WP_048165429.1",
        "ACS90033.1",
    ]
    mats = ProteinInfo.from_ncbi(mat_accessions)

    # Create pairwise alignments between all sequences
    alignments = Alignment.from_sequences(mats, aligner=PairwiseAligner)

    # Create a network
    network = SequenceNetwork(
        sequences=mats,
        pairwise_alignments=alignments,
        weight="identity",
        threshold=0.9,
        dimensions=3,
        color="taxonomy_id",
    )

    # Visualize the network
    network.visualize()
    ```

## Network Analysis

Upon the `SequenceNetwork` is instantiated the `graph` property is created. This property is a `networkx` graph object that can be used to perform network analysis. For example, the `degree()` method can be used to calculate the degree of each node in the network.
