import asyncio
import json
from logging import Logger
from typing import List

import nest_asyncio
from pyeed.adapter.primary_db_adapter import AsyncParamRequester, AsyncRequester
from pyeed.dbconnect import DatabaseConnector
from pyeed.fetch.dbsort import DBPattern, SortIDs
from pyeed.fetch.ncbiproteinmapper import NCBIProteinMapper
from pyeed.fetch.pdbmapper import PDBMapper
from pyeed.fetch.taxonomymapper import TaxonomyMapper
from rich.console import Console
from rich.progress import Progress


class ProteinFetcher:
    def __init__(self, ids: List[str], db: DatabaseConnector):
        self.ids = ids
        self.db = db
        nest_asyncio.apply()

    async def fetch(self, **console_kwargs):
        """
        Fetches protein data from various databases based on the provided IDs.

        Parameters:
            force_terminal (bool): Whether to force the use of a terminal
            for progress tracking.

        Returns:
            List[ProteinRecord]: A list of ProteinRecord objects containing the fetched protein data.

        Raises:
            Exception: If there is an error during the fetching process.
        """

        db_entries = SortIDs.sort(self.ids)
        param_requester = None

        with Progress(
            console=Console(**console_kwargs),
        ) as progress:
            requesters = []
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
                            n_concurrent=5,
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

                elif db_name == DBPattern.PDB.name:
                    task_id = progress.add_task(
                        f"Requesting sequences from {db_name}...", total=len(db_ids)
                    )

                    query_string = """
                    query {
                      entry(entry_id: "SEQUENCE_ID") {
                        polymer_entities {
                          rcsb_id
                          rcsb_polymer_entity_container_identifiers {
                            reference_sequence_identifiers {
                              database_accession
                              database_name
                            }
                          }
                          rcsb_entity_source_organism {
                            ncbi_taxonomy_id
                          }
                          entity_poly {
                            pdbx_seq_one_letter_code
                          }
                          rcsb_polymer_entity_feature {
                            type
                            feature_id
                            feature_positions {
                              beg_seq_id
                              end_seq_id
                            }
                          }
                          polymer_entity_instances {
                            rcsb_id
                            rcsb_polymer_instance_feature {
                              name
                              feature_positions {
                                beg_seq_id
                                end_seq_id
                              }
                            }
                          }
                        }
                      }
                    }
                    """

                    query = {"query": query_string}

                    param_requester = AsyncParamRequester(
                        ids=db_ids,
                        url="https://data.rcsb.org/graphql",
                        params=query,
                        task_id=task_id,
                        progress=progress,
                        batch_size=1,
                        rate_limit=50,
                        n_concurrent=20,
                    )

            responses = await asyncio.gather(
                *[requester.make_request() for requester in requesters]
            )

            if param_requester:
                pdb_entries = await param_requester.make_request()
                pdb_entries = [PDBMapper().map_pdb_data(entry) for entry in pdb_entries]
                pdb_entries = [item for sublist in pdb_entries for item in sublist]

            # map data to objects
            ncbi_responses, uniprot_response = self.identify_data_source(responses)

            ncbi_entries = NCBIProteinMapper().map(ncbi_responses)

            uniprot_entries = [
                UniprotMapper().map(*resp) for resp in uniprot_response.values()
            ]

            uniprot_entries.extend(ncbi_entries)
            if param_requester:
                uniprot_entries.extend(pdb_entries)

            # get taxonomy data
            unique_tax_ids = []
            for entry in uniprot_entries:
                if entry.organism and entry.organism.taxonomy_id not in unique_tax_ids:
                    unique_tax_ids.append(entry.organism.taxonomy_id)

            task_id = progress.add_task(
                "Requesting taxonomy data from EBI...", total=len(unique_tax_ids)
            )

            tax_requester = AsyncRequester(
                ids=unique_tax_ids,
                url="https://www.ebi.ac.uk/ena/browser/api/xml/",
                task_id=task_id,
                progress=progress,
                batch_size=1,
                rate_limit=50,
                n_concurrent=20,
            )

            taxonomies = await tax_requester.make_request()

            taxonomies = [
                tax for tax in taxonomies if "status" not in tax and "404" not in tax
            ]

            progress.update(task_id, completed=len(unique_tax_ids))

            # map taxonomy data to objects
            organisms = [TaxonomyMapper().map(entry) for entry in taxonomies]

            for entry in uniprot_entries:
                for organism in organisms:
                    if not entry.organism:
                        continue
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
        ncbi = None
        for response in responses:
            if response[0].startswith("LOCUS"):
                ncbi = response
            elif response[0].startswith("{"):
                uniprot[DBPattern.INTERPRO.name] = [
                    json.loads(entry) for entry in response
                ]
            elif response[0].startswith("["):
                uniprot[DBPattern.UNIPROT.name] = [
                    json.loads(entry)[0] for entry in response
                ]
            else:
                Logger.warning(f"Response could not be mapped to mapper: {response[0]}")

        if not ncbi:
            ncbi = []

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


if __name__ == "__main__":
    import asyncio

    from rich.progress import Progress

    ids = ["6VXX", "7NHM", "5L2G"]

    res = asyncio.run(ProteinFetcher(ids).fetch())
