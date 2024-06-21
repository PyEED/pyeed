import neo4j


class DBConnect:
    """Class to interact with neo4j database."""

    def __init__(
        self,
        url: str,
        db_name: str,
        user: str | None = None,
        password: str | None = None,
    ):
        self.url = url
        self.db_name = db_name
        self.user = user
        self.password = password

    def _get_driver(self):
        return neo4j.GraphDatabase.driver(
            self.url,
            auth=(self.user, self.password),
        )

    def add_node(self, node):
        return None

    def add_edge(self, node1, node2):
        return None
