import asyncio

import nest_asyncio

from pyeed.dbconnect import DatabaseConnector
from pyeed.fetch.mapper import UniprotToPyeed
from pyeed.fetch.requester import PrimaryDBRequester


class Pyeed:
    def __init__(
        self,
        uri: str,
        user: str | None = None,
        password: str | None = None,
    ):
        self.db = DatabaseConnector(uri, user, password)

    def fetch_from_primary_db(self, ids: list[str]):
        """
        Fetches sequences and corresponding annotations from primary sequence databases
        and adds them to local database.
        """
        nest_asyncio.apply()

        if isinstance(ids, str):
            ids = [ids]

        params_template = {
            "format": "json",
        }

        requester = PrimaryDBRequester(
            ids=ids,
            ids_attr_name="accession",
            url="https://www.ebi.ac.uk/proteins/api/proteins",
            rate_limit=10,
            n_concurrent=5,
            batch_size=1,
            data_mapper=UniprotToPyeed(),
            progress=None,
            task_id=None,
            request_params=params_template,
        )

        asyncio.run(requester.make_request())


if __name__ == "__main__":
    eedb = Pyeed("bolt://127.0.0.1:7687")
    eedb.db._wipe_database()

    eedb.fetch_from_primary_db(["P12345", "P67890", "P05062"])
    print(eedb.db.stats())
