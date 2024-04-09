# Clustering Sequences

`ProteinInfo` objects can be clustered using the `MMSeqs2` class. [MMSeqs2](https://github.com/soedinglab/mmseqs2) (Many-against-Many sequence searching) is a software suite to search and cluster huge protein and nucleotide sequence sets. PyEED provides a simple interface to the `easy-cluster` method of MMSeqs2. Thereby, the minimum sequence identity and the minimum coverage can be set to cluster the sequences. The `easy-cluster()` method returns an instance of ``MMSeqs2``. The created clusters can be accessed using the `clusters` attribute of the `MMSeqs2` object. Each cluster has a `representative` attribute that contains the representative sequence of the cluster and a `members` attribute that contains all the sequences in the cluster.

!!! example "Clustering Sequences"
    ```py
    from pyeed.core import ProteinInfo
    from pyeed.containers import MMSeqs2

    # Accessions of different methionine adenyltransferases
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
    mats = ProteinInfo.get_ids(mat_accessions)

    # Cluster the sequences
    clusterer = MMSeqs2.easy_cluster(
        sequences=mats,
        coverage=0.9,
        identity=0.6,
    )
    ```

The `easy-cluster()` method returns an instance of ``MMSeqs2``. The created clusters can be accessed using the `clusters` attribute of the `MMSeqs2` object. Each cluster has a `representative` attribute that contains the representative sequence of the cluster and a `members` attribute that contains all the sequences in the cluster. The clusters are sorted by descending member count.

!!! example "Accessing the clusters"
    ```py
    # Get the first cluster
    cluster = clusterer.clusters[0]

    # Get the representative sequence of the cluster
    representative = cluster.representative

    # Get the members of the cluster
    members = cluster.members
    ```