# this is the first general test for the fetching a protein and loaidng it into the database
import logging

from pyeed import Pyeed
from pyeed.model import GOAnnotation, Protein

LOGGER = logging.getLogger(__name__)


class TestProteinFetech:

    # this function is run before each test
    def set_up_test(self):
        LOGGER.info("Setting up test")

        uri = "bolt://127.0.0.1:7687"
        user = "neo4j"
        password = "12345678"

        # Create a Pyeed object, automatically connecting to the database
        eedb = Pyeed(uri, user, password)

        # For testing purposes, we will wipe the database and remove all constraints
        eedb.db._wipe_database()
        eedb.db._remove_db_constraints(user, password)

        # DB connector is an attribute of the Pyeed object, type `DatabaseConnector`
        print(eedb.db)
        print(type(eedb.db))

        # The first time the pyeed database is initialized, we need to create the constraints which are defined in the pyeed graph model
        eedb.db._initialize_db_constraints(user=user, password=password)

    def test_uniprot_fetch(self):
        LOGGER.info("Running test_uniprot_fetch")
        self.set_up_test()

        # Create a Pyeed object
    
