from pyeed.core import ProteinRecord
from pyeed.network import SequenceNetwork


class TestNetworkGraphBuild:
    def test_general_build_networkx(self):
        # check if it can be created and does the basic job
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
        mats = ProteinRecord.get_ids(mat_accessions)
        # Create a network
        network = SequenceNetwork(
            sequences=mats,
            weight="identity",
            dimensions=2,
        )

    def test_graph_build(self):
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
        assert len(list(network.network.nodes)) == len(mats)
