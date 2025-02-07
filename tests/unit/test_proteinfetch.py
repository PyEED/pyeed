# # this is the first general test for the fetching a protein and loading it into the database
# import logging

# from pyeed.main import Pyeed

# logging.basicConfig(
#     level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
# )
# LOGGER = logging.getLogger(__name__)


# class TestProteinFetech:
#     # this function is run before each test
#     def set_up_test(self):
#         LOGGER.info("Setting up test")

#         uri = "bolt://127.0.0.1:7687"
#         user = "neo4j"
#         password = "12345678"

#         # Create a Pyeed object, automatically connecting to the database
#         eedb = Pyeed(uri, user, password)

#         # For testing purposes, we will wipe the database and remove all constraints
#         eedb.db.wipe_database()
#         eedb.db.remove_db_constraints(user=user, password=password)

#         # DB connector is an attribute of the Pyeed object, type `DatabaseConnector`
#         LOGGER.info(f"Database stats: {eedb.db.stats()}")

#         # The first time the pyeed database is initialized, we need to create the constraints which are defined in the pyeed graph model
#         eedb.db.initialize_db_constraints(user=user, password=password)

#         self.eedb = eedb

#     def test_uniprot_fetch(self):
#         LOGGER.info("Running test_uniprot_fetch")
#         self.set_up_test()

#         # ids
#         ids = [
#             "P04182",
#             "Q6QDP7",
#         ]

#         # Fetch proteins from primary database
#         self.eedb.fetch_from_primary_db(ids)

#         # Check that the proteins were fetched
#         LOGGER.info(f"Database stats: {self.eedb.db.stats()}")
#         LOGGER.info("Test complete")

#         ids_in_db = self.eedb.getProteins()

#         LOGGER.info(f"Proteins in database: {ids_in_db}")

#         assert len(ids_in_db) == 2
#         # check if all ids are in the database
#         for id in ids:
#             assert id in [protein["accession_id"] for protein in ids_in_db]
#         # check if all ids in the database are in the ids
#         for protein in ids_in_db:
#             assert protein["accession_id"] in ids

#     def test_ncbi_fetch(self):
#         LOGGER.info("Running test_ncbi_fetch")
#         self.set_up_test()

#         # ids
#         ids = [
#             "AAL29438.1",
#             "HBQ2613975.1",
#             "EKW4005960.1",
#             "EJG7116187.1",
#             "AMM70781.1",
#             "HCO3480053.1",
#             "HAI5030310.1",
#             "WP_000027057.1",
#             "WP_215748091.1",
#             "WP_261627585.1",
#         ]

#         # Fetch proteins from primary database
#         self.eedb.fetch_from_primary_db(ids, db="NCBI")

#         # Check that the proteins were fetched
#         LOGGER.info(f"Database stats: {self.eedb.db.stats()}")
#         LOGGER.info("Test complete")

#         ids_in_db = self.eedb.getProteins()

#         LOGGER.info(f"Proteins in database: {ids_in_db}")

#         assert len(ids_in_db) == len(ids)
#         # check if all ids are in the database
#         for id in ids:
#             assert id in [protein["accession_id"] for protein in ids_in_db]
#         # check if all ids in the database are in the ids
#         for protein in ids_in_db:
#             assert protein["accession_id"] in ids

#     def test_ncbi_fetch_and_dna(self):
#         LOGGER.info("Running test_ncbi_fetch")
#         self.set_up_test()

#         # ids
#         ids = [
#             "AAL29438.1",
#             "HBQ2613975.1",
#             "EKW4005960.1",
#             "EJG7116187.1",
#             "AMM70781.1",
#             "HCO3480053.1",
#             "HAI5030310.1",
#             "WP_000027057.1",
#             "WP_215748091.1",
#             "WP_261627585.1",
#         ]

#         # Fetch proteins from primary database
#         self.eedb.fetch_from_primary_db(ids, db="NCBI")

#         # Check that the proteins were fetched
#         LOGGER.info(f"Database stats: {self.eedb.db.stats()}")
#         LOGGER.info("Test complete")

#         ids_in_db = self.eedb.getProteins()

#         LOGGER.info(f"Proteins in database: {ids_in_db}")

#         self.eedb.fetchRemoteCodingSequences()

#         assert len(ids_in_db) == len(ids)
#         # check if all ids are in the database
#         for id in ids:
#             assert id in [protein["accession_id"] for protein in ids_in_db]
#         # check if all ids in the database are in the ids
#         for protein in ids_in_db:
#             assert protein["accession_id"] in ids
