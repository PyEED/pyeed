# Creating Sequence Networks

A `SequenceNetwork` object can be created using a list of ProteinRecord objects. Those sequences are then used to create a alignment, based on this alignment each ProteinRecord objects represents a node. The edges can than be created based on a weight, e.g. 'identity', but custom weights can be introduced.
With the threshold mode and the threshold value, a threshold is then set to for example hide all edges with an identity score below 0.8.
The final network can be visualized and also loaded into cytoscape for further settings. Moreover it also can be used in maplotlib to plot, if intrested in custom styles.

## Visualization

``` py
    mat_accessions = [
        "MBP1912539.1",
        "SEV92896.1",
        "MBO8174569.1",
        "WP_042680787.1",
    ]
    mats = ProteinRecord.get_ids(mat_accessions)
    # Create network
    network = SequenceNetwork(
        sequences=mats,
        weight="identity",
    )

    network.create_graph()
    network.visualize()
```

Exporting the network in cytoscape is done the following way:
``` py
    import py4cytoscape as p4c

    # transfer the network to cytoscape
    netowork.create_cytoscape_graph(
        threshold=0.75,
        column_name="class",
    )

    # plot the network
    p4c.notebook_export_show_image()
``` 

The networkx object on which SequenceNetwork is based on can be extracted using the following command.

``` py
    graph_network = network.network
```

This then could be used to plot using matplotlib.

``` py

    import networkx as nx
    import matplotlib.pyplot as plt

    # Plotting of the network
    pos = nx.spring_layout(graph_network, weight='identity', iterations=100, seed=18)
    plt.figure(figsize=(19,9))
    nx.draw_networkx(graph_network, pos=pos, with_labels=True, node_color=c, node_size=s,
                    font_color='Black',font_size='6', font_weight='bold', edge_color='grey', alpha=0.5, width=1)
    plt.axis('off')
    plt.show()
``` 

## Network Analysis

Upon the `SequenceNetwork` is instantiated the `graph` property is created. This property is a `networkx` graph object that can be used to perform network analysis. For example, the `degree()` method can be used to calculate the degree of each node in the network.
