# Clustering Sequences

`ProteinInfo` objects can be clustered using the `MMSeqs2` class. [MMSeqs2](https://github.com/soedinglab/mmseqs2) (Many-against-Many sequence searching) is a software suite to search and cluster huge protein and nucleotide sequence sets. PyEED provides a simple interface to the `easy-cluster` method of MMSeqs2. Thereby, the minimum sequence identity and the minimum coverage can be set to cluster the sequences. The `easy-cluster()` method returns a list of `ProteinInfo` objects that represent the individual clusters with the specified identity and coverage.

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
    representatives = MMSeqs2.easy_cluster(
        sequences=mats,
        coverage=0.9,
        identity=0.6,
    )
    ```