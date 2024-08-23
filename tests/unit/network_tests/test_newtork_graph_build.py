import json
import pytest
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

    def test_graph_build(self):
        mat_accessions = [
            "MBP1912539.1",
            "SEV92896.1",
            "MBO8174569.1",
            "WP_042680787.1",
        ]
        mats = ProteinRecord.get_ids(mat_accessions)
        # Create network
        network = SequenceNetwork(sequences=mats)

        assert len(list(network.network.nodes)) == len(mats)

    def test_graph_save(self):
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

        # this is the data as a string, here one could write this to a file
        data_string = network.model_dump_json()
        # to read in the data again
        data_json = json.loads(data_string)
        # to load it into a Sequence Network object
        network2 = SequenceNetwork(**data_json)
                
        assert len(data_json['network']['nodes']) == len(mats)
        assert network2.network.number_of_nodes() == network.network.number_of_nodes()
        assert network2.network.number_of_edges() == network.network.number_of_edges()

    def test_graph_network_input_interger(self):
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

        # this is the data as a string, here one could write this to a file
        data_string = network.model_dump_json()
        # to read in the data again
        data_json = json.loads(data_string)

        data_json['network'] = 1

        # expecte an error
        with pytest.raises(ValueError):
            network2 = SequenceNetwork(**data_json)

