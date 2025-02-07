from abc import ABC, abstractmethod
from typing import Dict, Generic, List, NamedTuple, Optional, TypeVar

import aiometer
import tenacity
from httpx import AsyncClient, Limits, RequestError, Response, TimeoutException
from loguru import logger
from rich.progress import Progress, TaskID

# Define a type variable bound to our mapper interface
T = TypeVar("T", bound="PrimaryDBMapper")


class PrimaryDBMapper(ABC):
    """
    Abstract base class for mapping data from a primary sequence database
    to the pyeed graph object model and saving it to the database.
    """

    @abstractmethod
    def add_to_db(self, response: Response) -> None:
        """
        Maps the data from the primary database response to the pyeed model
        and saves it to the database.

        Args:
            response (Response): The HTTP response from the primary database.
        """
        ...


class RequestPayload(NamedTuple):
    """Container for HTTP request payload information."""

    client: AsyncClient
    url: str
    params: Dict[str, str]


class PrimaryDBAdapter(Generic[T]):
    """
    Orchestrates asynchronous HTTP GET requests to a primary sequence database.
    It uses an injected mapper (data_mapper) to process responses and save data.
    """

    def __init__(
        self,
        ids: List[str],
        id_param_name: str,
        url: str,
        rate_limit: int,
        max_concurrent: int,
        batch_size: int,
        data_mapper: T,
        timeout: int = 120,
        progress: Optional[Progress] = None,
        task_id: Optional[TaskID] = None,
        request_params: Optional[Dict[str, str]] = None,
    ):
        self.original_ids = ids  # the original list of IDs
        self.id_param_name = id_param_name
        self.url = url
        self.batch_size = batch_size
        self.rate_limit = rate_limit
        self.max_concurrent = max_concurrent
        self.data_mapper = data_mapper
        self.timeout = timeout
        self.request_params = request_params.copy() if request_params else {}

        # Create batches if batch_size > 1
        self.ids = self._create_batches() if batch_size > 1 else ids

        # Setup progress bar if not provided
        if progress is None:
            self.progress = Progress(disable=True)
            self.task_id = self.progress.add_task(
                "Requesting data...", total=len(self.ids)
            )
        else:
            self.progress = progress
            if task_id is not None:
                self.task_id = task_id

    def _create_batches(self) -> List[str]:
        """
        Groups IDs into comma-separated batch strings.
        Returns:
            List[str]: Batches of IDs as comma-separated strings.
        """
        return [
            ",".join(self.original_ids[i : i + self.batch_size])
            for i in range(0, len(self.original_ids), self.batch_size)
        ]

    def _create_payload(self, client: AsyncClient, id_value: str) -> RequestPayload:
        """
        Creates a RequestPayload using the provided client and an ID value.
        Args:
            client (AsyncClient): The HTTP client.
            id_value (str): The ID to include in the request parameters.
        Returns:
            RequestPayload: A named tuple with client, URL, and parameters.
        """
        params = self.request_params.copy()
        params[self.id_param_name] = id_value
        return RequestPayload(client=client, url=self.url, params=params)

    @tenacity.retry(
        wait=tenacity.wait_fixed(1),
        stop=tenacity.stop_after_attempt(3),
        retry=tenacity.retry_if_exception_type((RequestError, TimeoutException)),
    )
    async def _fetch_response(self, payload: RequestPayload) -> Response:
        """
        Sends an asynchronous HTTP GET request using the provided payload.
        Args:
            payload (RequestPayload): The payload for the request.
        Returns:
            Response: The HTTP response.
        """
        logger.debug(
            f"Sending request to {payload.url} with parameters: {payload.params}"
        )
        return await payload.client.get(
            payload.url, params=payload.params, timeout=self.timeout
        )

    def _update_progress(self) -> None:
        """
        Advances the progress bar by one unit.
        """
        if self.progress and self.task_id is not None:
            self.progress.update(self.task_id, advance=1)

    async def execute_requests(self) -> None:
        """Executes the asynchronous HTTP GET requests concurrently."""
        async with AsyncClient(
            limits=Limits(max_connections=self.max_concurrent)
        ) as client:
            logger.info(f"Starting requests for {len(self.ids)} batches.")
            payloads = [self._create_payload(client, id_value) for id_value in self.ids]
            logger.debug(f"Prepared {len(payloads)} request payloads.")

            async with aiometer.amap(
                self._fetch_response,
                payloads,
                max_per_second=self.rate_limit,
                max_at_once=self.max_concurrent,
            ) as response_coroutines:
                async for response in response_coroutines:
                    if response.status_code == 200 and response.content:
                        self.data_mapper.add_to_db(response)
                    else:
                        logger.error(
                            f"Request failed with status {response.status_code}: {response.text}"
                        )
                    self._update_progress()
