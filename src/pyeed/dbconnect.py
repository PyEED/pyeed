import os
import subprocess

from neo4j import Driver, GraphDatabase
from neomodel import db as neomodel_db


class DatabaseConnector:
    def __init__(self, uri: str, user: str | None, password: str | None):
        """
        Initializes the connection to the Neo4j database using a self-managed driver.
        """
        self._uri = uri
        self.driver = self._get_driver(uri, user, password)
        neomodel_db.set_connection(driver=self.driver)  # patch db for neomodel

        if not self.constraints_exist():
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

    def execute_read(self, query: str, parameters=None) -> list[dict]:
        """
        Executes a read (MATCH) query using the Neo4j driver.

        Args:
            query (str): The Cypher query to execute.
            parameters (dict): A dictionary of parameters to pass to the query.

        Returns:
            list[dict]: The result of the query as a list of dictionaries.
        """
        with self.driver.session() as session:
            return session.execute_read(self._run_query, query, parameters)

    def execute_write(self, query: str, parameters=None):
        """
        Executes a write (CREATE, DELETE, etc.) query using the Neo4j driver.

        Args:
            query (str): The Cypher query to execute.
            parameters (dict): A dictionary of parameters to pass to the query.
        """
        with self.driver.session() as session:
            return session.execute_write(self._run_query, query, parameters)

    def stats(self) -> dict:
        """
        Returns the number of nodes and relationships in the database.

        Returns:
            dict: The number of nodes and relationships in the database.
        """
        node_count_query = "MATCH (n) RETURN count(n) AS node_count"
        relationship_count_query = (
            "MATCH ()-[r]->() RETURN count(r) AS relationship_count"
        )

        node_count = self.execute_read(node_count_query)[0]["node_count"]
        relationship_count = self.execute_read(relationship_count_query)[0][
            "relationship_count"
        ]

        return {"nodes": node_count, "relationships": relationship_count}

    def initialize_db_constraints(
        self,
        user: str | None,
        password: str | None,
        models_path: str = "model.py",
    ):
        """
        Run the neomodel_install_labels script to set up indexes and constraints on labels
        of Object-Graph Mapping (OGM) models.

        Args:
            user (str): The username for the Neo4j database.
            password (str): The password for the Neo4j database.
            models_path (str): The path to the models file. Defaults to "model.py".
        """
        # set the path to the models file
        # work from the path of this file
        models_path = os.path.join(os.path.dirname(__file__), models_path)

        try:
            # Construct connection string based on whether user/password are provided
            if user and password and self._uri.startswith("bolt"):
                connection_url = f"bolt://{user}:{password}@{self._uri.split('//')[1]}"
                print("the connection url is", connection_url)
            elif user and password and self._uri.startswith("neo4j+s"):
                connection_url = (
                    f"neo4j+s://{user}:{password}@{self._uri.split('//')[1]}"
                )
                print("the connection url is", connection_url)
            else:
                connection_url = self._insert_after_second_slash(
                    self._uri, "neo4j:neo4j@"
                )
                print("the connection url is", connection_url)

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

    def constraints_exist(self) -> bool:
        """
        Check and if constraints exist in the database. Return True if constraints exist.

        Returns:
            bool: True if constraints exist in the database, False otherwise.
        """
        query = """
        SHOW CONSTRAINTS YIELD name, type
        RETURN count(*) AS constraint_count
        """

        results = self.execute_read(query)
        return True if results[0]["constraint_count"] > 0 else False

    def remove_db_constraints(
        self,
        user: str | None,
        password: str | None,
    ):
        """
        Run the neomodel_remove_labels script to drop all indexes and constraints
        from labels in the Neo4j database.

        Args:
            user (str): The username for the Neo4j database.
            password (str): The password for the Neo4j database
        """
        try:
            if user and password:
                connection_url = f"bolt://{user}:{password}@{self._uri.split('//')[1]}"
            else:
                connection_url = self._insert_after_second_slash(
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

    def generate_model_diagram(
        self,
        models_path: str = "pyeed/model.py",
    ):
        """Generates a arrows json file representing the model diagram.

        Args:
            models_path (str, optional): The path to the models file. Defaults to "pyeed/model.py".
        """
        subprocess.run(
            [
                "neomodel_generate_diagram",
                models_path,
            ]
        )

    def wipe_database(self):
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

    @property
    def node_properties(self) -> list[dict]:
        """
        Returns a list of dictionaries containing the node labels and their properties.
        """
        node_properties_query = """
            CALL apoc.meta.data()
            YIELD label, other, elementType, type, property
            WHERE NOT type = "RELATIONSHIP" AND elementType = "node"
            WITH label AS nodeLabels, collect(property) AS properties
            RETURN {labels: nodeLabels, properties: properties} AS output
            """

        return self.execute_read(node_properties_query)

    @property
    def relationship_properties(self) -> list[dict]:
        """
        Returns a list of dictionaries containing the relationship types and their properties.
        """
        rel_properties_query = """
            CALL apoc.meta.data()
            YIELD label, other, elementType, type, property
            WHERE NOT type = "RELATIONSHIP" AND elementType = "relationship"
            WITH label AS nodeLabels, collect(property) AS properties
            RETURN {type: nodeLabels, properties: properties} AS output
            """

        return self.execute_read(rel_properties_query)

    @property
    def relationships(self) -> list[dict]:
        """
        Returns a list of dictionaries containing the source node label, relationship type, and target node label.
        """
        rel_query = """
            CALL apoc.meta.data()
            YIELD label, other, elementType, type, property
            WHERE type = "RELATIONSHIP" AND elementType = "node"
            RETURN {source: label, relationship: property, target: other} AS output
            """

        return self.execute_read(rel_query)

    @staticmethod
    def _insert_after_second_slash(uri: str, to_insert: str) -> str:
        """Inserts a string after the second '//' in a URI.

        Args:
            uri (str): The URI to modify.
            to_insert (str): The string to insert.

        Returns:
            str: The modified URI.
        """
        scheme, rest = uri.split("//", 1)
        return f"{scheme}//{to_insert}{rest}"
