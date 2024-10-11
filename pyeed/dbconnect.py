import subprocess

from neo4j import Driver, GraphDatabase
from neomodel import db as neomodel_db

from pyeed.model import Protein


class DatabaseConnector:
    def __init__(self, uri: str, user: str | None, password: str | None):
        """
        Initializes the connection to the Neo4j database using a self-managed driver.
        """
        self._uri = uri
        self.driver = self._get_driver(uri, user, password)
        neomodel_db.set_connection(driver=self.driver)  # patch db for neomodel

        if not self._constraints_exist():
            print(
                "Pyeed Graph Object Mapping constraints not defined. Use _install_labels() to set up model constraints."
            )
        print("📡 Connected to database.")

    def close(self):
        """
        Closes the connection to the Neo4j database.
        """
        self.driver.close()
        print("🔌 Connection closed.")

    def execute_read(self, query: str, parameters=None):
        """
        Executes a read (MATCH) query using the Neo4j driver directly.
        """
        with self.driver.session() as session:
            return session.execute_read(self._run_query, query, parameters)

    def execute_write(self, query: str, parameters=None):
        """
        Executes a write (CREATE, DELETE, etc.) query using the Neo4j driver directly.
        """
        with self.driver.session() as session:
            return session.execute_write(self._run_query, query, parameters)

    def add_protein(self, protein_record: Protein):
        """
        Placeholder for adding a Protein to the database via Neomodel.
        """
        # Here you can add logic to store protein_record using Neomodel models
        pass

    def _initialize_db_constraints(
        self,
        user: str | None,
        password: str | None,
        models_path: str = "pyeed/model.py",
    ):
        """
        Run the neomodel_install_labels script to set up indexes and constraints on labels
        of Object-Graph Mapping (OGM) models.
        """
        try:
            # Construct connection string based on whether user/password are provided
            if user and password:
                connection_url = f"bolt://{user}:{password}@{self._uri.split('//')[1]}"
            else:
                connection_url = self.insert_after_second_slash(
                    self._uri, "neo4j:neo4j@"
                )

            subprocess.run(
                [
                    "neomodel_install_labels",
                    models_path,
                    "--db",
                    connection_url,
                ],
                check=True,
            )
            print(
                "✅ Databse constraints and indexes set up according to Pyeed Graph Object Model."
            )
        except subprocess.CalledProcessError as e:
            print(f"Failed to install labels: {str(e)}")

    def _constraints_exist(self) -> bool:
        """Check and if constraints exist in the database. Return True if constraints exist."""
        query = """
        SHOW CONSTRAINTS YIELD name, type
        RETURN count(*) AS constraint_count
        """

        results = self.execute_read(query)
        return True if results[0]["constraint_count"] > 0 else False

    def _remove_db_constraints(
        self,
        user: str | None,
        password: str | None,
    ):
        """
        Run the neomodel_remove_labels script to drop all indexes and constraints
        from labels in the Neo4j database.
        """
        try:
            if user and password:
                connection_url = f"bolt://{user}:{password}@{self._uri.split('//')[1]}"
            else:
                connection_url = self.insert_after_second_slash(
                    self._uri, "neo4j:neo4j@"
                )

            subprocess.run(
                [
                    "neomodel_remove_labels",
                    "--db",
                    connection_url,
                ],
                check=True,
            )
            print("All constraints and indexes have been removed from the database.")
        except subprocess.CalledProcessError as e:
            print(f"Failed to remove labels: {str(e)}")

    def _wipe_database(self):
        """
        Deletes all nodes and relationships in the database.
        """
        delete_query = """
        MATCH (n)
        DETACH DELETE n
        """
        self.execute_write(delete_query)
        print("All data has been wiped from the database.")

    @staticmethod
    def _run_query(tx, query, parameters):
        """
        Executes a Cypher query in the provided transaction.
        """
        result = tx.run(query, parameters)
        return [record.data() for record in result]

    @staticmethod
    def _get_driver(uri: str, user: str | None, password: str | None) -> Driver:
        """
        Creates a new Neo4j driver instance.
        """
        auth = (user, password) if user and password else None
        return GraphDatabase.driver(uri, auth=auth)

    @staticmethod
    def insert_after_second_slash(uri: str, to_insert: str) -> str:
        # Split the string at '//' into two parts
        scheme, rest = uri.split("//", 1)

        # Insert the new content after the second '//'
        return f"{scheme}//{to_insert}{rest}"
