import io
from abc import ABC, abstractmethod
from typing import Any, Coroutine, Generic, NamedTuple, TypeVar

import aiometer
import tenacity
from Bio import SeqIO
from httpx import (
    AsyncClient,
    Limits,
    RequestError,
    Response,
    TimeoutException,
)
from loguru import logger
from rich.progress import Progress, TaskID

T = TypeVar("T", bound="PrimaryDBtoPyeed")


class PrimaryDBtoPyeed(ABC):
    """
    Abstract base class for mapping data from a primary sequence database to the pyeed
    graph object model and saving it to the database.
    """

    @abstractmethod
    def add_to_db(self, data: Any) -> None:
        """Abstract method for mapping data from a primary sequence database to the pyeed
        graph object model and saving it to the database.

        Args:
            data (Any): The data to be mapped and saved to the database.
        """
        pass


class RequestPayload(NamedTuple):
    """Holds the request client, URL, and parameters for an HTTP GET request."""

    client: AsyncClient
    url: str
    params: dict[str, str]


class PrimaryDBAdapter(Generic[T]):
    """
    Orchestrates the asynchronous HTTP GET requests to a primary sequence database.
    Mapper classes are injected to map the responses to the pyeed graph object model and
    save them to the database.
    """

    def __init__(
        self,
        ids: list[str],
        ids_attr_name: str,
        url: str,
        rate_limit: int,
        n_concurrent: int,
        batch_size: int,
        data_mapper: T,
        timeout: int = 120,
        progress: Progress | None = None,
        task_id: TaskID | None = None,
        request_params: dict[str, str] = {},
    ):
        self.ids = ids
        self.ids_attr_name = ids_attr_name
        self.url = url
        self.batch_size = batch_size
        self.rate_limit = rate_limit
        self.n_concurrent = n_concurrent
        self.progress = progress
        self.task_id = task_id
        self.data_mapper = data_mapper
        self.timeout = timeout
        self.request_params = request_params

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
        Groups the IDs into batches of the specified batch size.

        Returns:
            list[str]: The list of batches, where each batch is a comma-separated
            string of IDs.
        """
        batches = []
        for i in range(0, len(self.ids), self.batch_size):
            batch = self.ids[i : i + self.batch_size]
            batch_string = ",".join(map(str, batch))
            batches.append(batch_string)
        return batches

    def build_request_payload(self, client: AsyncClient, id_: str) -> RequestPayload:
        """Combines the client, URL, and parameters into a RequestPayload object.
        Adds the id with the key specified by ids_attr_name to the request parameters.

        Args:
            client (AsyncClient): AsyncClient object for making HTTP requests
            id_ (str): ID to be added to the request parameters

        Returns:
            RequestPayload: RequestPayload object with the client, URL, and parameters
        """
        params = self.request_params.copy()
        params[self.ids_attr_name] = id_

        return RequestPayload(client, self.url, params=params)

    @tenacity.retry(
        wait=tenacity.wait_fixed(1),
        stop=tenacity.stop_after_attempt(3),
        retry=tenacity.retry_if_exception_type((RequestError, TimeoutException)),
    )
    async def send_request(
        self,
        args: RequestPayload,
    ) -> Coroutine[None, None, Response]:
        """
        Sends an asynchronous HTTP GET request to the specified URL using the provided
        AsyncClient.
        """
        client = args.client
        url = args.url
        params = args.params

        logger.debug(f"Sending request to {url} with parameters: {params}")
        return client.get(url, params=params, timeout=self.timeout)

    async def make_request(self):
        """
        Makes asynchronous HTTP GET requests to the specified URL using the provided
        AsyncClient, handling rate limiting and concurrency.
        """

        def update_progress():
            if self.progress and self.task_id:
                self.progress.update(self.task_id, advance=1)  # type: ignore

        async with AsyncClient(
            limits=Limits(max_connections=self.n_concurrent),
        ) as client:
            logger.info(f"Making requests with ids list: {self.ids}")
            # Build the list of request arguments (this prepares the coroutine tasks)
            requests = [self.build_request_payload(client, id) for id in self.ids]

            logger.debug(
                f"Sending {len(self.ids)} requests in batches of {self.batch_size}"
            )

            # Using aiometer to handle rate-limiting and concurrency
            async with aiometer.amap(
                self.send_request,
                requests,
                max_per_second=self.rate_limit,
                max_at_once=self.n_concurrent,
            ) as response_coroutines:
                async for response_coroutine in response_coroutines:
                    res = await response_coroutine
                    sanitized_response = self.sanitize_response(res)
                    logger.debug(f"Received response: {sanitized_response}")
                    [self.map_and_add_to_db(entry) for entry in sanitized_response]

                    update_progress()

    def sanitize_response(self, response: Response) -> list[dict[str, Any]]:
        """
        Sanitizes the response from the HTTP GET request by checking the status code
        and formatting the JSON response as a list of dictionaries.

        Returns:
            Optional[List[Dict[str, Any]]]: The JSON response as a list of dictionaries,
            or None if the response is invalid.
        """
        if response.status_code != 200:
            logger.warning(
                f"Request to {response.url} failed with status code {response.status_code}"
            )
            return []

        try:
            logger.debug(f"Response content: {type(response.content)}")

            # here we need to identify from where the response is coming from and parse it accordingly
            if response.content.startswith(b"LOCUS"):
                return SeqIO.parse(io.StringIO(response.content.decode()), "gb")
            elif response.content.startswith(b"{"):
                None
            elif response.content.startswith(b"["):
                response_json = response.json()
                if not response_json:
                    logger.warning(f"Empty response from {response.url}")
                    return []

                # If the response is a dictionary, wrap it in a list
                if isinstance(response_json, dict):
                    response_json = [response_json]

                # Ensure the response is a list of dictionaries
                if not isinstance(response_json, list) or not all(
                    isinstance(item, dict) for item in response_json
                ):
                    logger.warning(f"Unexpected response format from {response.url}")
                    return []

                return response_json

            else:
                logger.warning(f"Response could not be mapped to mapper: {response}")

        except ValueError as e:
            logger.warning(f"Failed to parse JSON response from {response.url}: {e}")
            return []

        return response_json

    def map_and_add_to_db(self, response: dict[str, Any] | None):
        """
        Handles the response from the HTTP GET request by passing it to the data mapper.
        This adds the mapped data to the database.
        """

        if response is None:
            return None
        self.data_mapper.add_to_db(response)
