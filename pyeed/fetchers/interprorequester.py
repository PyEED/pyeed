import asyncio
from typing import List
import aiohttp


class InterProRequester:

    def __init__(self, ids: List[str], rate_limit_per_second: int = 30):
        self.ids = ids
        self.rate_limit_per_second = rate_limit_per_second

    async def make_request(self):
        async with aiohttp.ClientSession() as session:
            tasks = [self.fetch(session, id) for id in self.ids]
            return await asyncio.gather(*tasks)

    async def fetch(self, session, id):
        url = (
            "https://www.ebi.ac.uk/interpro/api/entry/all/protein/UniProt/"
            + id
            + "?format=json"
        )
        async with session.get(url) as response:
            data = await response.json()
            await asyncio.sleep(1 / self.rate_limit_per_second)
            return data
