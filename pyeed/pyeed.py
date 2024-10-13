import asyncio

import nest_asyncio

from pyeed.adapter.primary_db_adapter import PrimaryDBAdapter
from pyeed.adapter.uniprot_mapper import UniprotToPyeed
from pyeed.dbconnect import DatabaseConnector


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

        adapter = PrimaryDBAdapter(
            ids=ids,
            ids_attr_name="accession",
            url="https://www.ebi.ac.uk/proteins/api/proteins",
            rate_limit=10,
            n_concurrent=5,
            batch_size=5,
            data_mapper=UniprotToPyeed(),
            progress=None,
            task_id=None,
            request_params=params_template,
        )

        asyncio.run(adapter.make_request())


if __name__ == "__main__":
    eedb = Pyeed("bolt://127.0.0.1:7687")

    search = False
    if search:
        eedb.db._wipe_database()

        eedb.fetch_from_primary_db(
            [
                "P04182",
                "Q6QDP7",
                "P04182",
                "P29758",
                "A0A851UXD9",
                "A0A8C6HVU6",
                "A0A8C6GQ10",
                "A0A1U7QEB0",
                "A0A6I9L5L6",
                "G3HVE0",
                "A0A8J6G992",
                "A0A8C6W4W5",
                "A0A8B9YUY7",
                "L8I4V3",
                "A0A6P3IYQ1",
                "A0A452EKJ3",
                "A0A6P5B7Q0",
                "F1MYG0",
                "A0A5J5MK22",
                "A0A6J0Y425",
                "Q3ZCF5",
            ]
        )

    print(eedb.db.stats())
