import asyncio
import logging
from typing import List

import nest_asyncio
from pyeed.adapter.primary_db_adapter import AsyncRequester
from rich.console import Console
from rich.progress import Progress

from .ncbidnamapper import NCBIDNAMapper

LOGGER = logging.getLogger(__name__)


class DNAFetcher:
    def __init__(self, ids: List[str]):
        self.ids = ids
        nest_asyncio.apply()

    async def fetch(self, **console_kwargs):
        """
        Fetches DNA data from various databases based on the provided IDs.

        Parameters:
            force_terminal (bool): Whether to force the use of a terminal
            for progress tracking.

        Returns:
            List[dnarecord]: A list of dnaRecord objects containing the fetched dna data.

        Raises:
            Exception: If there is an error during the fetching process.

        """
        # right now in the first batch version we just fetch from NCBI
        param_requester = None

        with Progress(
            console=Console(**console_kwargs),
        ) as progress:
            requesters: List[AsyncRequester] = []

            #
            task_id = progress.add_task(
                "Requesting sequences from NCBI...", total=len(self.ids)
            )
            requesters.append(
                AsyncRequester(
                    ids=self.ids,
                    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&retmode=text&rettype=genbank&id=",
                    task_id=task_id,
                    progress=progress,
                    batch_size=10,
                    rate_limit=2,
                    n_concurrent=5,
                )
            )

            responses = await asyncio.gather(
                *[requester.make_request() for requester in requesters]
            )

            # in case of multiple databases, identify the source of the data
            ncbi_responses, uniprot_response = self.identify_data_source(responses)

            # map data to objects
            ncbi_entries = NCBIDNAMapper().map(responses=ncbi_responses)

            return ncbi_entries

    def identify_data_source(self, responses: List[str]) -> tuple:
        """
        Identifies the source of the data based on the response content.
        """
        ncbi_responses = []
        uniprot_response = []

        for response in responses:
            if response[0].startswith("LOCUS"):
                ncbi_responses.append(response)
            else:
                uniprot_response.append(response)

        return ncbi_responses, uniprot_response
