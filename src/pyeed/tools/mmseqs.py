import logging
from typing import List, Optional

import httpx
from loguru import logger
from pyeed.dbconnect import DatabaseConnector
from pyeed.tools.abstract_tool import AbstractTool, ServiceURL
from pyeed.tools.datamodels.mmseqs_model import Cluster


class MMSeqs(AbstractTool):
    """
    Class for MMSeqs search with MMSeqs2 Tool
    """

    def __init__(self):
        super().__init__()
        self._service_url = ServiceURL.MMSEQS.value

    def extract_output_data(self, result) -> List[Cluster]:
        """
        Extracts the output data from the ClustalOmega container.

        Returns:
            MultiSequenceAlignment: The alignment result.
        """

        cleaned_text = result.text.encode().decode(
            "unicode_escape"
        )  # Decode Unicode escapes
        cleaned_text = cleaned_text.replace(
            ">", ""
        ).splitlines()  # Remove '>' and split into lines

        if cleaned_text:  # Ensure the list is not empty
            cleaned_text.pop(-1)  # Remove the last element if needed

        # Use list comprehension to filter elements
        cleaned_text = [element for element in cleaned_text if len(element) < 25]

        print(cleaned_text)

        # get all the representatives and the belonging represented sequences

        clusters = []
        seq_counter = 0

        while seq_counter < len(cleaned_text) - 1:  # Iterate over all sequences
            # If the current and next sequence are the same, it's a representative sequence
            if cleaned_text[seq_counter] == cleaned_text[seq_counter + 1]:
                representative_id = cleaned_text[seq_counter]
                represented_ids = [representative_id]  # Include itself as represented

                # Add the cluster
                clusters.append(
                    Cluster(
                        representative_id=representative_id,
                        represented_ids=represented_ids,
                    )
                )

                seq_counter += 2
            else:
                # Skip sequences that do not follow the doubling pattern
                seq_counter += 1

        return clusters

    def run_service(self, data, min_seq_id, coverage, cov_mode) -> httpx.Response:
        """_summary_

        Args:
            query (str): The query is either the sequence as it would be entered in a FASTA file or a id of a sequence in the database.
            result_filename (str): The name of the file to store the results.
            min_seq_id (float): min sequence identity for which a cluster is created
            coverage (float): min coverage for which a cluster is created
            cov_mode (int): coverage mode- choose between 0, 1, 2

        Returns:
            httpx.Response: The response object containing the MMSeqs results.
        """
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)

        if cov_mode not in [0, 1, 2]:
            raise ValueError("cov_mode must be 0, 1, or 2")

        json_request = {
            "query": data,
            "min_seq_id": min_seq_id,
            "coverage": coverage,
            "cov_mode": cov_mode,
        }

        self._service_url = ServiceURL.MMSEQS.value

        try:
            # Log before making the request
            logger.info(
                f"Sending POST request to {self._service_url} with payload: {json_request}"
            )
            return httpx.post(self._service_url, json=json_request, timeout=6000)

        except httpx.ConnectError as connect_error:
            context = connect_error.__context__
            if context and hasattr(context, "args"):
                error_number = context.args[0].errno
                if error_number == 8 or error_number == -3:
                    self._service_url = "http://localhost:8001/easycluster"
                    try:
                        return httpx.post(
                            self._service_url, json=json_request, timeout=6000
                        )
                    except httpx.ConnectError:
                        raise httpx.ConnectError(
                            "PyEED Docker Service not running"
                        ) from None
            print(connect_error)
            raise httpx.ConnectError("PyEED Docker Service not running") from None

    def easycluster(
        self,
        query: str,
        min_seq_id: float = 0.5,
        coverage: float = 0.8,
        cov_mode: int = 0,
        dbConnector: Optional[DatabaseConnector] = None,
    ) -> List[Cluster]:
        """
        Clusters the sequences in the fasta file with the provided parameters
        Args:
            query (str): The query is either the sequence as it would be entered in a FASTA file or a id of a sequence in the database.
            min_seq_id (float): min sequence identity for which a cluster is created
            coverage (float): min coverage for which a cluster is created
            cov_mode (int): coverage mode- choose between 0, 1, 2

        Returns:
            List[Cluster]: A list of Cluster objects containing the representative and all sequences.
        """

        if dbConnector:
            query_run = """
            MATCH (p:Protein)
            WHERE p.accession_id IN $ids
            RETURN p.accession_id AS accession_id, p.sequence AS sequence
            """

            # Execute the query for all IDs in the `query` list
            proteins = dbConnector.execute_read(query_run, {"ids": query})

            # Check if any proteins were found
            if not proteins:
                raise ValueError("No proteins found with the provided IDs.")

            # Convert the results into FASTA format
            fasta_sequences = []
            for protein in proteins:
                id = protein["accession_id"]
                protein_seq = protein["sequence"]
                fasta = f">{id}\n{protein_seq}"
                fasta_sequences.append(fasta)

            # Combine all FASTA sequences into a single string if needed
            fasta_output = "\n".join(fasta_sequences)
            query = fasta_output

        result = self.run_service(query, min_seq_id, coverage, cov_mode)
        logger.info(f"MMSeqs2 search completed with status code: {result.status_code}")

        if result.status_code != 200:
            logger.error(f"Error: {result.json()}")

        return self.extract_output_data(result)
