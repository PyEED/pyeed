import io
from typing import Optional

import httpx
import pandas as pd
from loguru import logger
from pyeed.dbconnect import DatabaseConnector
from pyeed.tools.abstract_tool import AbstractTool, ServiceURL


class Blast(AbstractTool):
    """
    Class for Blast search with biopython
    """

    def __init__(self):
        super().__init__()
        self._service_url = ServiceURL.BLASTP.value

    def run_service(
        self, query, db, evalue, outfmt, num_threads, max_target_seqs, protein=True
    ):
        """
        Run the BLAST service with the provided parameters.

        Args:
            query (str): The query is either the sequence as it would be entered in a FASTA file or a id of a sequence in the database.
            db (str): The database to search against.
            evalue (float): E-value threshold for BLAST hits.
            outfmt (int): Output format for BLAST results.
            num_threads (int): Number of threads to use for the BLAST search.
            max_target_seqs (int): Maximum number of target sequences to retrieve.
            protein (bool): If True, perform a BLASTP search. If False, perform a BLASTN search.

        Returns:
            httpx.Response: The response object containing the BLAST results.

        """

        json_request = {
            "query": query,
            "db": db,
            "evalue": evalue,
            "outfmt": outfmt,
            "num_threads": num_threads,
            "max_target_seqs": max_target_seqs,
        }

        if protein:
            self._service_url = ServiceURL.BLASTP.value
        else:
            self._service_url = ServiceURL.BLASTN.value

        try:
            return httpx.post(self._service_url, json=json_request, timeout=6000)

        except httpx.ConnectError as connect_error:
            context = connect_error.__context__
            if context and hasattr(context, "args"):
                error_number = context.args[0].errno
                if error_number == 8 or error_number == -3:
                    if protein:
                        self._service_url = "http://localhost:6001/blastp"
                    else:
                        self._service_url = "http://localhost:6001/blastn"
                    try:
                        return httpx.post(
                            self._service_url, json=json_request, timeout=6000
                        )
                    except httpx.ConnectError:
                        raise httpx.ConnectError(
                            "PyEED Docker Service is not running."
                        ) from None

            print(connect_error)
            raise httpx.ConnectError("PyEED Docker Service is not running.") from None

    def blastp(
        self,
        query,
        db,
        evalue=0.001,
        outfmt=10,
        num_threads=50,
        max_target_seqs=50,
        dbConnector: Optional[DatabaseConnector] = None,
    ) -> pd.DataFrame:
        """
        Perform a BLASTP search using the provided parameters.

        Args:
            query (str): The query is either the sequence as it would be entered in a FASTA file or a id of a sequence in the database.
            db (str): The database to search against.
            evalue (float): E-value threshold for BLAST hits. Default is 0.001.
            outfmt (int): Output format for BLAST results. Default is 10 (CSV).
            num_threads (int): Number of threads to use for the BLAST search. Default is 50.
            max_target_seqs (int): Maximum number of target sequences to retrieve. Default is 50.
            dbConnector (DatabaseConnector): A DatabaseConnector object to connect to the database. Default is None.

        Returns:
            pd.DataFrame: A pandas DataFrame containing the BLAST results.
        """
        # check if the query is an id or a sequence
        if dbConnector:
            query_run = """
            MATCH (p:Protein)
            WHERE p.accession_id IN $ids
            RETURN p.accession_id AS accession_id, p.sequence AS sequence
            """
            protein = dbConnector.execute_read(query_run, {"ids": [query]})

            if len(protein) == 0:
                raise ValueError("No protein found with the provided id.")

            protein_seq = protein[0]["sequence"]
            id = protein[0]["accession_id"]

            # convert in fasta format
            query = f">{id}\n{protein_seq}"

        # now we will run the run service
        result = self.run_service(
            query, db, evalue, outfmt, num_threads, max_target_seqs
        )
        logger.info(f"BLASTP search completed with status code {result.status_code}")

        if result.status_code != 200:
            logger.error(f"Error: {result.json()}")

        # create a dataframe from the result
        df = pd.read_csv(io.StringIO(result.json()["result"]))
        # add the rigth header
        header = "Query ID,Subject ID,% Identity,Alignment Length,Mismatches,Gap Opens,Query Start,Query End,Subject Start,Subject End,E-value,Bit Score".split(
            ","
        )
        df.columns = header

        return df

    def blastn(
        self,
        query,
        db,
        evalue=0.001,
        outfmt=10,
        num_threads=50,
        max_target_seqs=50,
        dbConnector: Optional[DatabaseConnector] = None,
    ) -> pd.DataFrame:
        """
        Perform a BLASTN search using the provided parameters.

        Args:
            query (str): The query is either the sequence as it would be entered in a FASTA file or a id of a sequence in the database.
            db (str): The database to search against.
            evalue (float): E-value threshold for BLAST hits. Default is 0.001.
            outfmt (int): Output format for BLAST results. Default is 10 (CSV).
            num_threads (int): Number of threads to use for the BLAST search. Default is 50.
            max_target_seqs (int): Maximum number of target sequences to retrieve. Default is 50.
            dbConnector (DatabaseConnector): A DatabaseConnector object to connect to the database. Default is None.

        Returns:
            pd.DataFrame: A pandas DataFrame containing the BLAST results.
        """
        # check if the query is an id or a sequence
        if dbConnector:
            query_run = """
            MATCH (n:DNA)
            WHERE n.accession_id IN $ids
            RETURN n.accession_id AS accession_id, n.sequence AS sequence
            """
            dnas = dbConnector.execute_read(query_run, {"ids": [query]})

            if len(dnas) == 0:
                raise ValueError("No DNA found with the provided id.")

            dna_seq = dnas[0]["sequence"]
            id = dnas[0]["accession_id"]

            # convert in fasta format
            query = f">{id}\n{dna_seq}"

        # now we will run the run service
        result = self.run_service(
            query, db, evalue, outfmt, num_threads, max_target_seqs, protein=False
        )
        logger.info(f"BLASTN search completed with status code {result.status_code}")

        if result.status_code != 200:
            logger.error(f"Error: {result.json()}")

        # create a dataframe from the result
        df = pd.read_csv(io.StringIO(result.json()["result"]))
        # add the rigth header
        header = "Query ID,Subject ID,% Identity,Alignment Length,Mismatches,Gap Opens,Query Start,Query End,Subject Start,Subject End,E-value,Bit Score".split(
            ","
        )
        df.columns = header

        return df
