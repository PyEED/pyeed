import json
import pytest

from pyeed.network import SequenceNetwork
from pyeed.core import ProteinInfo

class TestNetworkDataframe:

    def test_general(self):
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
        mats = ProteinInfo.get_ids(mat_accessions)
        # Create a network
        network = SequenceNetwork(
            sequences=mats,
            weight="identity",
            threshold=0.9,
            dimensions=2,
            color="taxonomy_id",
        )

