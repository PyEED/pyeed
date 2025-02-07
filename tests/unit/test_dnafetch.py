# # this is the first general test for the fetching a dna and loading it into the database
# import logging

# from pyeed import Pyeed

# logging.basicConfig(
#     level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
# )
# LOGGER = logging.getLogger(__name__)


# class TestDNAFetech:
#     # this function is run before each test
#     def set_up_test(self):
#         LOGGER.info("Setting up test")

#         uri = "bolt://127.0.0.1:7687"
#         user = "neo4j"
#         password = "12345678"

#         # Create a Pyeed object, automatically connecting to the database
#         eedb = Pyeed(uri, user, password)

#         # For testing purposes, we will wipe the database and remove all constraints
#         eedb.db._wipe_database()
#         eedb.db._remove_db_constraints(user, password)

#         # DB connector is an attribute of the Pyeed object, type `DatabaseConnector`
#         LOGGER.info(f"Database stats: {eedb.db.stats()}")

#         # The first time the pyeed database is initialized, we need to create the constraints which are defined in the pyeed graph model
#         eedb.db._initialize_db_constraints(user=user, password=password)

#         self.eedb = eedb

#     def test_fetch_dna_from_ncbi(self):
#         LOGGER.info("Running dna fetch test")
#         self.set_up_test()

#         ids = [
#             "AF397067.1",
#             "AF397068.1",
#         ]

#         # Fetch dna from primary database
#         self.eedb.fetch_nucleotide_from_db(ids)

#         # Check that the proteins were fetched
#         LOGGER.info(f"Database stats: {self.eedb.db.stats()}")
#         LOGGER.info("Test complete")

#         ids_in_db = self.eedb.getDNAs()

#         LOGGER.info(f"DNAS in database: {ids_in_db}")

#         assert len(ids_in_db) == 2
#         # check if all ids are in the database
#         for id in ids:
#             assert id in [protein["accession_id"] for protein in ids_in_db]
#         # check if all ids in the database are in the ids
#         for protein in ids_in_db:
#             assert protein["accession_id"] in ids
