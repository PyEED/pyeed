import pandas as pd
import networkx as nx

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
        )
        network.create_graph()

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

        network.create_graph()

        assert len(list(network.network.nodes)) == len(mats)

    def test_graph_naming_dic(self):
        mat_accessions = [
            "MBP1912539.1",
            "SEV92896.1",
            "MBO8174569.1",
            "WP_042680787.1",
        ]

        naming_dic = {
            "MBP1912539": "MAT1",
            "SEV92896": "MAT2",
            "MBO8174569": "MAT3",
            "WP_042680787": "MAT4",
        }

        mats = ProteinRecord.get_ids(mat_accessions)
        # Create network
        network = SequenceNetwork(
            sequences=mats,
            weight="identity",
            naming_dict=naming_dic,
        )

        network.create_graph(naming_dict=naming_dic)

        # convert networkx to pandas
        data = dict(network.network.nodes(data=True))

        assert len(list(network.network.nodes)) == len(mats)
        assert 'MAT1' in list(data.keys())