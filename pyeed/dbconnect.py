from neo4j import Driver, GraphDatabase

from pyeed.model import ProteinRecord


class DatabaseConnector:
    def __init__(self, uri: str, user: str | None, password: str | None):
        self.driver = self._get_driver(uri, user, password)
        print("📡 Connected to database.")

    def close(self):
        self.driver.close()

    def execute_read(self, query, parameters=None):
        with self.driver.session() as session:
            return session.execute_read(self._run_query, query, parameters)

    def execute_write(self, query, parameters=None):
        with self.driver.session() as session:
            return session.execute_write(self._run_query, query, parameters)

    def add_protein(self, protein_record: ProteinRecord):
        pass

    def _wipe_database(self):
        delete_query = """
        MATCH (n)
        DETACH DELETE n
        """
        self.execute_write(delete_query)
        print("All data has been wiped from the database.")

    @staticmethod
    def _run_query(tx, query, parameters):
        result = tx.run(query, parameters)
        return [record.data() for record in result]

    @staticmethod
    def _get_driver(uri: str, user: str | None, password: str | None) -> Driver:
        auth = (user, password) if user and password else None
        return GraphDatabase.driver(uri, auth=auth)


# Example Usage
if __name__ == "__main__":
    db = DatabaseConnector("bolt://127.0.0.1:7687", None, None)
    try:
        # Sample nested JSON
        json_data = {
            "label": "ProteinRecord",
            "id": "P12345",
            "name": "Protein A",
            "sequence": "MSEQ1234",
            "organism": {
                "label": "Organism",
                "taxonomy_id": 9606,
                "name": "Homo sapiens",
            },
            "structure": [
                {
                    "label": "Structure",
                    "id": "STR123",
                    "description": "Protein structure 1",
                },
                {
                    "label": "Structure",
                    "id": "STR456",
                    "description": "Protein structure 2",
                },
            ],
        }

        # Create nodes and relationships based on the JSON
        db.create_nodes_from_json(json_data)

        # Verify: Get the number of nodes in the database
        count_query = "MATCH (n) RETURN count(n) AS node_count"
        node_count_result = db.execute_read(count_query)
        node_count = node_count_result[0]["node_count"]
        print(f"Number of nodes in the database: {node_count}")

    finally:
        db.close()
