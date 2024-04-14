import json
import asyncio
import logging
import nest_asyncio
from typing import List
from rich.console import Console
from rich.progress import Progress

from pyeed.fetch.dbsort import SortIDs, DBPattern
from pyeed.fetch.uniprotmapper import UniprotMapper
from pyeed.fetch.ncbiproteinmapper import NCBIProteinMapper
from pyeed.fetch.taxonomymapper import TaxonomyMapper
from pyeed.fetch.requester import AsyncRequester


LOGGER = logging.getLogger(__name__)


class ProteinFetcher:

    def __init__(self, ids: List[str]):
        self.ids = ids
        # self.ncbi_key = ncbi_key #TODO: Add NCBI key to NCBI requester
        nest_asyncio.apply()

    async def fetch(self, force_terminal: bool = False):
        """
        Fetches protein data from various databases based on the provided IDs.

        Parameters:
            force_terminal (bool): Whether to force the use of a terminal
            for progress tracking.

        Returns:
            List[ProteinInfo]: A list of ProteinInfo objects containing the fetched protein data.

        Raises:
            Exception: If there is an error during the fetching process.

        """
        db_entries = SortIDs.sort(self.ids)

        with Progress(console=Console(force_terminal=force_terminal)) as progress:
            requesters: List[AsyncRequester] = []
            for db_name, db_ids in db_entries.items():

                if db_name == DBPattern.UNIPROT.name:
                    task_id = progress.add_task(
                        f"Requesting sequences from {db_name}...", total=len(db_ids)
                    )
                    requesters.append(
                        AsyncRequester(
                            ids=db_ids,
                            url="https://www.ebi.ac.uk/proteins/api/proteins?format=json&accession=",
                            task_id=task_id,
                            progress=progress,
                            batch_size=1,
                            rate_limit=5,
                            n_concurrent=20,
                        )
                    )

                    task_id = progress.add_task(
                        "Requesting domains from INTERPRO...", total=len(db_ids)
                    )
                    requesters.append(
                        AsyncRequester(
                            ids=db_ids,
                            url="https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/",
                            task_id=task_id,
                            progress=progress,
                            batch_size=1,
                            rate_limit=10,
                            n_concurrent=10,
                        )
                    )

                elif db_name == DBPattern.NCBI.name:
                    task_id = progress.add_task(
                        f"Requesting sequences from {db_name}...", total=len(db_ids)
                    )
                    requesters.append(
                        AsyncRequester(
                            ids=db_ids,
                            url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=text&rettype=genbank&id=",
                            task_id=task_id,
                            progress=progress,
                            batch_size=10,
                            rate_limit=2,
                            n_concurrent=20,
                        )
                    )

                elif db_name == DBPattern.INTERPRO.name:
                    task_id = progress.add_task(
                        f"Requesting domains from {db_name}...", total=len(db_ids)
                    )
                    requesters.append(
                        AsyncRequester(
                            ids=db_ids,
                            url="https://www.ebi.ac.uk/interpro/api/entry/interpro/",
                            task_id=task_id,
                            progress=progress,
                            rate_limit=10,
                            n_concurrent=10,
                            batch_size=1,
                        )
                    )

            responses = await asyncio.gather(
                *[requester.make_request() for requester in requesters]
            )

            # map data to objects
            ncbi_response, uniprot_response = self.identify_data_source(responses)

            ncbi_entries = NCBIProteinMapper().map(ncbi_response)
            uniprot_entries = [
                UniprotMapper().map(*resp) for resp in uniprot_response.values()
            ]

            uniprot_entries.extend(ncbi_entries)

            # get taxonomy data
            unique_tax_ids = list(
                set([entry.organism.taxonomy_id for entry in uniprot_entries])
            )

            task_id = progress.add_task(
                "Requesting taxonomy data from EBI...", total=len(unique_tax_ids)
            )

            tax_requester = AsyncRequester(
                ids=unique_tax_ids,
                url="https://www.ebi.ac.uk/ena/browser/api/xml/",
                task_id=task_id,
                progress=progress,
                batch_size=1,
                rate_limit=10,
                n_concurrent=20,
            )

            taxonomies = await tax_requester.make_request()

            # map taxonomy data to objects
            organisms = [TaxonomyMapper().map(entry) for entry in taxonomies]

            for entry in uniprot_entries:
                for organism in organisms:
                    if entry.organism.taxonomy_id == organism.taxonomy_id:
                        entry.organism = organism

            return uniprot_entries

    def identify_data_source(self, responses: List[List[str]]):
        """
        Identifies the data source for each response in the list of responses.

        Parameters:
            responses (List[List[str]]): A list of responses, where each response
            is a list of strings.

        Returns:
            Tuple[List[str], dict]: A tuple containing the NCBI response and
            a dictionary of UniProt responses sorted by database pattern.

        """
        uniprot = {}
        for response in responses:
            if response[0].startswith("LOCUS"):
                ncbi = response
                print("NCBI detected")
            elif response[0].startswith("{"):
                print("INTERPRO detected")
                uniprot[DBPattern.INTERPRO.name] = [
                    json.loads(entry) for entry in response
                ]
            elif response[0].startswith("["):
                print("UNIPROT detected")
                uniprot[DBPattern.UNIPROT.name] = [
                    json.loads(entry)[0] for entry in response
                ]
            else:
                LOGGER.warning(f"Response could not be mapped to mapper: {response[0]}")

        if uniprot:
            uniprot_dict = self.sort_uniprot_by_id(uniprot)
        else:
            uniprot_dict = {}

        return ncbi, uniprot_dict

    def sort_uniprot_by_id(self, uniprot: dict) -> dict:
        """
        Sorts the UniProt responses by their ID and returns a dictionary.

        Parameters:
            uniprot (dict): A dictionary containing UniProt responses sorted
            by database pattern.

        Returns:
            dict: A dictionary where the keys are the UniProt IDs and the
            values are the corresponding UniProt responses.

        """
        uniprot_dict = {}
        for uniprot_response in uniprot[DBPattern.UNIPROT.name]:
            target_id = uniprot_response["accession"].upper()
            uniprot_dict[target_id] = uniprot_response

            for interpro_response in uniprot[DBPattern.INTERPRO.name]:
                if (
                    target_id
                    == interpro_response["results"][0]["proteins"][0][
                        "accession"
                    ].upper()
                ):
                    uniprot_dict[target_id] = (uniprot_response, interpro_response)

        return uniprot_dict
