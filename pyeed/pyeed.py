import nest_asyncio

from pyeed.dbconnect import DatabaseConnector


class pyeed:
    def __init__(
        self,
        uri: str,
        user: str | None = None,
        password: str | None = None,
    ):
        self.db = DatabaseConnector(uri, user, password)
        _with_db: bool = True

    def fetch_from_primary_db(self, ids: list[str]):
        """
        Fetches sequences and corresponding annotations from primary sequence databases
        and adds them to local database.
        """
        nest_asyncio.apply()

        if isinstance(ids, str):
            ids = [ids]
