from typing import Any, Dict, Generic, List, NamedTuple, Optional, TypeVar

import aiometer
import tenacity
from httpx import AsyncClient, Limits, Response
from loguru import logger as LOGGER
from rich.progress import Progress, TaskID

from pyeed.dbconnect import DatabaseConnector
from pyeed.fetch.mapper import PrimaryDBResponseMapper

T = TypeVar("T")


class RequestArgs(NamedTuple):
    """Holds the arguments for an HTTP GET request."""

    client: AsyncClient
    url: str
    params: Optional[Dict[str, str]] = None


class PrimaryDBRequester(Generic[T]):
    def __init__(
        self,
        ids: List[str],
        url: str,
        rate_limit: int,
        n_concurrent: int,
        batch_size: int,
        data_mapper: "PrimaryDBResponseMapper[T]",
        db: DatabaseConnector,
        progress: Optional[Progress] = None,
        task_id: Optional[TaskID] = None,
        params_template: Optional[Dict[str, str]] = None,
        use_params: bool = False,
    ):
        self.ids = ids
        self.url = url
        self.batch_size = batch_size
        self.rate_limit = rate_limit
        self.n_concurrent = n_concurrent
        self.progress = progress
        self.task_id = task_id
        self.data_mapper = data_mapper
        self.database_connector = db
        self.params_template = params_template
        self.use_params = use_params

        if self.batch_size > 1:
            self.ids = self.make_batches()

        if not self.progress:
            self._create_progress()

    def _create_progress(self):
        """
        Creates a dummy progress bar for tracking the progress of the HTTP
        requests if not provided.
        """
        self.progress = Progress(disable=True)
        self.task_id = self.progress.add_task("Requesting data...", total=len(self.ids))

    def make_batches(self) -> list[str]:
        """
        Creates batches of IDs for making HTTP requests.

        Returns:
            List[str]: The list of batches, where each batch is a comma-separated
            string of IDs.
        """
        batches = []
        for i in range(0, len(self.ids), self.batch_size):
            batch = self.ids[i : i + self.batch_size]
            batch_string = ",".join(map(str, batch))
            batches.append(batch_string)
        return batches

    def build_request_args(self, client: AsyncClient, id: str) -> RequestArgs:
        if self.use_params and self.params_template:
            params = self.params_template.copy()
            for key, value in params.items():
                params[key] = value.replace("SEQUENCE_ID", id)
            return RequestArgs(client, self.url, params=params)
        else:
            url = f"{self.url}{id}"
            return RequestArgs(client, url)

    @tenacity.retry(
        wait=tenacity.wait_fixed(1),  # Wait 1 second between retries
        stop=tenacity.stop_after_attempt(1),  # Retry up to 3 times
    )
    async def send_request(self, args: RequestArgs) -> Optional[Dict[str, Any]]:
        """
        Sends an asynchronous HTTP GET request to the specified URL using the provided
        AsyncClient.
        """
        client = args.client
        url = args.url
        params = args.params

        LOGGER.debug(f"Sending request to {url}")
        response = await client.get(url, params=params, timeout=120)

        LOGGER.debug(f"Received response from {url}. Code: {response.status_code}")

        if response.status_code != 200:
            LOGGER.warning(
                f"Request to {url} failed with status code {response.status_code}"
            )
            return None  # Early return if the request failed

        try:
            response_json = response.json()

            # Check if the response is a list and has elements
            if isinstance(response_json, list) and response_json:
                response_dict = response_json[0]
            elif isinstance(response_json, dict):
                response_dict = response_json
            else:
                LOGGER.warning(f"Unexpected response format from {url}")
                return None  # Return None if response format is unexpected

        except ValueError as e:
            LOGGER.error(f"Failed to parse JSON response from {url}: {str(e)}")
            return None  # Return None if JSON parsing fails

        return response_dict

        # Use data mapper if provided

    async def make_request(self) -> List[T]:
        """
        Makes asynchronous HTTP GET requests to the specified URL using the provided
        AsyncClient.

        Returns:
            List[Any]: The list of transformed data from the responses.
        """
        all_responses: list[T] = []

        async def update_progress(response: Response):
            if self.progress and self.task_id:
                self.progress.update(self.task_id, advance=1)  # type: ignore

        async with AsyncClient(
            event_hooks={"response": [update_progress]},
            limits=Limits(max_connections=self.n_concurrent),
        ) as client:
            LOGGER.debug(f"Creating {len(self.ids)} tasks")

            tasks = [self.build_request_args(client, id) for id in self.ids]

            LOGGER.debug(f"Sending {len(self.ids)} requests")
            async with aiometer.amap(
                self.send_request,
                tasks,
                max_per_second=self.rate_limit,
                max_at_once=self.n_concurrent,
            ) as responses:
                async for res in responses:
                    if res is not None:
                        data = await self.data_mapper.transform(res)

                        all_responses.append(data)

        return all_responses


if __name__ == "__main__":
    import asyncio

    from pyeed.fetch.mapper import UniprotMapper

    requester = PrimaryDBRequester(
        ids=["P12345", "P67890", "P54321"],
        url="https://www.ebi.ac.uk/proteins/api/proteins?format=json&accession=",
        rate_limit=10,
        n_concurrent=5,
        batch_size=1,
        data_mapper=UniprotMapper(),
        db=None,
        progress=None,
        task_id=None,
        params_template=None,
        use_params=False,
    )

    responses = asyncio.run(requester.make_request())
    seq = responses[0]

    with open("seq.json", "w") as f:
        f.write(seq.model_dump_json())
