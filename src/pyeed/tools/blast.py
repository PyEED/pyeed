import io
from typing import Literal, Optional

import httpx
import pandas as pd
from loguru import logger
from pyeed.dbconnect import DatabaseConnector
from pyeed.tools.abstract_tool import AbstractTool


class Blast(AbstractTool):
    def run_service(
        self,
        mode: Literal["blastp", "blastn"],
        query: str,
        db_path: str,
        db_name: str,
        evalue: float,
        max_target_seqs: int,
        num_threads: int,
    ) -> httpx.Response:
        """Run the BLAST service with the provided parameters.

        Args:
            query: Query sequence
            db: Path to BLAST database
            evalue: E-value threshold
            num_threads: Number of threads to use
            max_target_seqs: Maximum number of target sequences
            mode: BLAST mode ("blastp" or "blastn")

        Returns:
            Response from the BLAST service
        """
        if mode not in ["blastp", "blastn"]:
            raise ValueError("Invalid mode. Must be 'blastp' or 'blastn'.")
        logger.debug(f"Running {mode} search")

        data = {
            "mode": mode,
            "sequence": query,
            "db_path": db_path,
            "db_name": db_name,
            "evalue": evalue,
            "max_target_seqs": max_target_seqs,
            "num_threads": num_threads,
        }

        try:
            return httpx.post(
                "http://localhost:6001/blast",
                json=data,  # Send as JSON body, not query params
                timeout=6000,
            )
        except httpx.ConnectError as e:
            logger.error(f"Connection error: {e}")
            raise httpx.ConnectError("PyEED Docker Service not running") from e

    def blastp(
        self,
        query,
        db,
        num_threads,
        evalue=0.001,
        outfmt=6,
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
            MATCH (p:DNA)
            WHERE p.accession_id IN $ids
            RETURN p.accession_id AS accession_id, p.sequence AS sequence
            """
            dnas = dbConnector.execute_read(query_run, {"ids": [query]})

            if len(dnas) == 0:
                raise ValueError("No protein found with the provided id.")

            dna_seq = dnas[0]["sequence"]
            id = dnas[0]["accession_id"]

            # convert in fasta format
            query = f">{id}\n{dna_seq}"

        # now we will run the run service
        result = self.run_service(
            query, db, evalue, outfmt, num_threads, max_target_seqs
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


if __name__ == "__main__":
    seq = "MSEQVAAVAKLRAKASEAAKEAKAREAAKKLAEAAKKAKAKEAAKRAEAKLAEKAKAAKRAEAKAAKEAKRAAAKRAEAKLAEKAKAAK"
    blast = Blast()
    res = blast.run_service(
        mode="blastp",
        query=seq,
        db_path="/usr/local/bin/data/test_db/",
        db_name="protein_db",
        evalue=10,
        max_target_seqs=10,
        num_threads=4,
    )

    print(res.json())
