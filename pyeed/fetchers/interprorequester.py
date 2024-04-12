import asyncio
from typing import List
import httpx


class InterProRequester:

    def __init__(self, ids: List[str], rate_limit_per_second: int = 30):
        self.ids = ids
        self.rate_limit_per_second = rate_limit_per_second

    async def make_request(self):
        """
        Make a request to the InterPro API for each ID in the list of IDs.

        Returns:
            A list of JSON responses from the InterPro API for each ID.

        Raises:
            aiohttp.ClientError: If there is an error making the request.
            asyncio.TimeoutError: If the request times out.
        """

        async with httpx.AsyncClient() as session:
            tasks = [self.fetch(session, id) for id in self.ids]
            return await asyncio.gather(*tasks)

    async def fetch(self, session: httpx.AsyncClient, id):
        """
        Fetches data from the InterPro API for a given ID.
        """
        url = f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/UniProt/{id}?format=json"
        try:
            response = await session.get(url)
            if response.status_code == 204:
                print(f"No InterPro entry found for ID: {id}")
                return None

            data = response.json()
            await asyncio.sleep(1 / self.rate_limit_per_second)

            return data

        except httpx.RequestError as e:
            print(f"Request error occurred: {e}")
        except httpx.HTTPStatusError as e:
            print(f"HTTP status error occurred: {e}")
        finally:
            await response.aclose()
