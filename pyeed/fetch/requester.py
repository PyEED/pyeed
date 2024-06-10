import asyncio
import logging
from typing import Dict, List, NamedTuple, Optional

import aiometer
import tenacity
from httpx import AsyncClient, Limits, Response
from rich.progress import Progress, TaskID

LOGGER = logging.getLogger(__name__)


class RequestArgs(NamedTuple):
    """Holds the arguments for an HTTP GET request."""

    client: AsyncClient
    url: str
    params: Optional[dict] = None


class AsyncRequester:
    def __init__(
        self,
        ids: List[str],
        url: str,
        batch_size: int,
        rate_limit: int,
        n_concurrent: int,
        progress: Optional[Progress] = None,
        task_id: Optional[TaskID] = None,
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

    @tenacity.retry(
        wait=tenacity.wait_fixed(1),
        stop=tenacity.stop_after_attempt(3),
    )
    async def send_request(self, args: RequestArgs) -> str:
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
        response = await client.get(url, timeout=120)

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

    @tenacity.retry(
        wait=tenacity.wait_fixed(0.5),
        stop=tenacity.stop_after_attempt(3),
    )
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
                self.progress.update(self.task_id, advance=self.batch_size)  # type: ignore

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

    def make_batches(self) -> List[str]:
        """
        Creates batches of IDs for making HTTP requests.

        Returns:
            List[str]: The list of batches, where each batch is a comma-separated
            string of IDs.

        """

        batches = []
        for i in range(0, len(self.ids), self.batch_size):
            batch = self.ids[i : i + self.batch_size]
            if len(batch) > 1:
                batch_string = ",".join(batch)
            else:
                batch_string = str(batch[0])
            batches.append(batch_string)
        self.ids = batches
        return batches


class AsyncParamRequester:
    """Updated Requester utilizing parameters as dict for the request"""

    def __init__(
        self,
        params: Dict[str, str],
        url: str,
        ids: List[str],
        rate_limit: int,
        n_concurrent: int,
        progress: Optional[Progress] = None,
        task_id: Optional[TaskID] = None,
        batch_size: int = 1,
    ):
        self.params = params
        self.url = url
        self.ids = ids
        self.batch_size = batch_size
        self.rate_limit = rate_limit
        self.n_concurrent = n_concurrent
        self.progress = progress
        self.task_id = task_id

        if not self.progress:
            self._create_progress()

    def _create_progress(self):
        """
        Creates a dummy progress bar for tracking the progress of the HTTP
        requests if not provided.
        """

        self.progress = Progress(disable=True)
        self.task_id = self.progress.add_task("Requesting data...", total=len(self.ids))

    async def send_request(self, args: RequestArgs) -> str:
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
        params = args.params

        LOGGER.debug(f"Sending request to {url}")
        response = await client.get(url, params=params, timeout=120)

        LOGGER.debug(f"Received response from {url}. Code: {response.status_code}")

        if response.status_code != 200:
            LOGGER.warning(
                f"Request to {url} failed with status code {response.status_code}"
            )
            LOGGER.warning(f"Response: {response.text}")

        if response.status_code == 429:
            LOGGER.warning("Rate limit exceeded. Waiting for 1 second...")
            await asyncio.sleep(0.5)
            return await self.send_request(args)

        return response.text

    async def make_request(self) -> List[str]:
        """Handles the asynchronous HTTP GET and configures rate limits and progress bar."""

        all_responses = []

        async def update_progress(response: Response):
            if self.progress:
                self.progress.update(self.task_id, advance=self.batch_size)  # type: ignore

        async with AsyncClient(
            event_hooks={"response": [update_progress]},
            limits=Limits(max_connections=self.n_concurrent, keepalive_expiry=30),
        ) as client:
            tasks = []
            for id in self.ids:
                params = self.params.copy()
                params["query"] = params["query"].replace("SEQUENCE_ID", str(id))
                tasks.append(RequestArgs(client, self.url, params))

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
