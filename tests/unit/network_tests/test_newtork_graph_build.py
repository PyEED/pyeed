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

        assert network.network is not None


    def test_cytoscope(self):
        mat_accessions = [
            "MBP1912539.1",
            "SEV92896.1",
            "MBO8174569.1",
            "WP_042680787.1",
        ]
        mats = ProteinRecord.get_ids(mat_accessions)
        # Create a network
        network = SequenceNetwork(
            sequences=mats,
            weight="identity",
            dimensions=2,
        )
        threshhold = 0.85
        # now create a cytoscape graph
        # careful if at this point cytoscope is not installed and running in background test will fail
        network.create_cytoscape_graph(collection="tests", title="test_cytoscape", threshold=threshhold)
        # check if the hidden edeges are above the threshold
        assert [item[1] for item in list(network._get_edges_visibilities().items())].count(False) == 4

    def test_cytoscape_degree(self):
        mat_accessions = [
            "MBP1912539.1",
            "SEV92896.1",
            "MBO8174569.1",
            "WP_042680787.1",
        ]
        mats = ProteinRecord.get_ids(mat_accessions)
        # Create a network
        network = SequenceNetwork(
            sequences=mats,
            weight="identity",
            dimensions=2,
        )
        threshhold = 0.85
        # now create a cytoscape graph
        network.create_cytoscape_graph(collection="tests", title="test_cytoscape", threshold=threshhold)
        # check if each degree is correct
        # it shoudl be at 1
        assert network.network.nodes["MBP1912539.1"]["degree_with_threshold_{}".format(threshhold)] == 1
        assert network.network.nodes["SEV92896.1"]["degree_with_threshold_{}".format(threshhold)] == 1
        assert network.network.nodes["MBO8174569.1"]["degree_with_threshold_{}".format(threshhold)] == 1
        assert network.network.nodes["WP_042680787.1"]["degree_with_threshold_{}".format(threshhold)] == 1

    
    def test_cytoscape_nodes_size(self):
        mat_accessions = [
            "MBP1912539.1",
            "SEV92896.1",
            "MBO8174569.1",
            "WP_042680787.1",
        ]
        mats = ProteinRecord.get_ids(mat_accessions)
        # Create a network
        network = SequenceNetwork(
            sequences=mats,
            weight="identity",
            dimensions=2,
        )
        threshhold = 0.85
        # now create a cytoscape graph
        network.create_cytoscape_graph(collection="tests", title="test_cytoscape", threshold=threshhold)
        network.set_nodes_size(column_name="degree_with_threshold_{}".format(threshhold), min_size=20, max_size=100)
        network.color_nodes(column_name='species')

        
