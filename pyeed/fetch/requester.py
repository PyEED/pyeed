import asyncio
import logging
import aiometer
from typing import List, NamedTuple
from rich.progress import Progress, TaskID
from httpx import AsyncClient, Limits, Response


LOGGER = logging.getLogger(__name__)


class RequestArgs(NamedTuple):
    """Holds the arguments for an HTTP GET request."""

    client: AsyncClient
    url: str


class AsyncRequester:

    def __init__(
        self,
        ids: List[str],
        url: str,
        batch_size: int = None,
        rate_limit: int = None,
        n_concurrent: int = None,
        progress: Progress = None,
        task_id: TaskID = None,
    ):
        self.ids = ids
        self.url = url
        self.batch_size = batch_size
        self.rate_limit = rate_limit
        self.n_concurrent = n_concurrent
        self.progress = progress
        self.task_id = task_id

        if self.batch_size:
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

    async def send_request(self, args: RequestArgs):
        """
        Sends an asynchronous HTTP GET request to the specified URL using the provided
        AsyncClient.

        Parameters:
            args (RequestArgs): The arguments for the request, including the client and
            the URL.

        Returns:
            str: The response text from the request.
        """

        client = args.client
        url = args.url

        LOGGER.debug(f"Sending request to {url}")
        response = await client.get(url, timeout=30)

        LOGGER.debug(f"Received response from {url}. Code: {response.status_code}")

        if response.status_code != 200:
            LOGGER.warning(
                f"Request to {url} failed with status code {response.status_code}"
            )
            LOGGER.warning(f"Response: {response.text}")

        if response.status_code == 429:
            LOGGER.warning("Rate limit exceeded. Waiting for 1 second...")
            await asyncio.sleep(1)
            return await self.send_request(args)

        return response.text

    async def make_request(self) -> List[str]:
        """
        Makes asynchronous HTTP GET requests to the specified URL using the provided
        AsyncClient.

        Returns:
            List[str]: The response texts from the requests.

        Notes:
            - If the response status code is not 200, a warning message is logged.
            - If the response status code is 429 (rate limit exceeded), the method waits
            for 1 second and then retries the request.
        """

        all_responses = []

        async def update_progress(response: Response):
            if self.progress:
                self.progress.update(self.task_id, advance=self.batch_size)

        async with AsyncClient(
            event_hooks={"response": [update_progress]},
            limits=Limits(max_connections=self.n_concurrent),
        ) as client:
            LOGGER.debug(f"Creating {len(self.ids)} tasks")

            tasks = [RequestArgs(client, f"{self.url}{id}") for id in self.ids]

            LOGGER.debug(f"Sending {len(self.ids)} requests")
            async with aiometer.amap(
                self.send_request,
                tasks,
                max_per_second=self.rate_limit,
                max_at_once=self.n_concurrent,
            ) as responses:
                async for res in responses:
                    all_responses.append(res)

        return all_responses

    def make_batches(self):
        """
        Creates batches of IDs for making HTTP requests.

        Returns:
            List[str]: The list of batches, where each batch is a comma-separated
            string of IDs.

        """

        batches = []
        for i in range(0, len(self.ids), self.batch_size):
            batch = self.ids[i : i + self.batch_size]
            batch_string = ",".join(batch)
            batches.append(batch_string)

        self.ids = batches
        return batches
