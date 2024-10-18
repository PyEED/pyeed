# this is the first general test for the fetching a protein and loaidng it into the database
import logging

from pyeed import Pyeed
from pyeed.model import GOAnnotation, Protein

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
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
        LOGGER.info(f"Database stats: {eedb.db.stats()}")

        # The first time the pyeed database is initialized, we need to create the constraints which are defined in the pyeed graph model
        eedb.db._initialize_db_constraints(user=user, password=password)

        self.eedb = eedb

    def test_uniprot_fetch(self):
        LOGGER.info("Running test_uniprot_fetch")
        self.set_up_test()

        # ids 
        ids = [
            "P04182",
            "Q6QDP7",
        ]

        # Fetch proteins from primary database
        self.eedb.fetch_from_primary_db(ids)

        # Check that the proteins were fetched
        LOGGER.info(f"Database stats: {self.eedb.db.stats()}")
        LOGGER.info("Test complete")

        ids_in_db = self.eedb.getProteins()

        LOGGER.info(f"Proteins in database: {ids_in_db}")

        assert len(ids_in_db) == 2
        assert ids_in_db[0]['accession_id'] == "P04182"
        assert ids_in_db[1]['accession_id'] == "Q6QDP7"

    def test_ncbi_fetch(self):
        LOGGER.info("Running test_ncbi_fetch")
        self.set_up_test()

        # ids 
        ids = [
            "AAL29438",
        ]

        # Fetch proteins from primary database
        self.eedb.fetch_from_primary_db(ids, db="NCBI")

        # Check that the proteins were fetched
        LOGGER.info(f"Database stats: {self.eedb.db.stats()}")
        LOGGER.info("Test complete")

        ids_in_db = self.eedb.getProteins()

        LOGGER.info(f"Proteins in database: {ids_in_db}")

        assert len(ids_in_db) == 2
        assert ids_in_db[0]['accession_id'] == "NP_001356"
        assert ids_in_db[1]['accession_id'] == "NP_001356.1"


        

    
