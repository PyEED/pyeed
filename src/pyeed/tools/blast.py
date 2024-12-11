import os
import httpx
import subprocess
from Bio.Blast.Record import Blast as BlastRecord
from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastnCommandline
from pyeed.tools.abstract_tool import AbstractTool, ServiceURL
from Bio.Blast import NCBIXML
import io
import pandas as pd

class Blast(AbstractTool):
    """
        Class for Blast search with biopython
    """

    def __init__(self):
        super().__init__()
        self._service_url = ServiceURL.BLASTP.value

    def run_service(self, query, db, evalue, outfmt, num_threads, max_target_seqs):
        
        
        json_request = {
            'query': query,
            'db': db,
            'evalue': evalue,
            'outfmt': outfmt,
            'num_threads': num_threads,
            'max_target_seqs': max_target_seqs
        }
        
    
        try:
            return httpx.post(self._service_url, json=json_request, timeout=6000)

        except httpx.ConnectError as connect_error:
            context = connect_error.__context__
            if context and hasattr(context, "args"):
                error_number = context.args[0].errno
                if error_number == 8 or error_number == -3:
                    self._service_url = 'http://localhost:6001/blastp'
                    try:
                        return httpx.post(self._service_url, json=json_request, timeout=6000)
                    except httpx.ConnectError:
                        raise httpx.ConnectError(
                            "PyEED Docker Service is not running."
                        ) from None

            print(connect_error)
            raise httpx.ConnectError("PyEED Docker Service is not running.") from None

    BLAST_COLUMNS = [
    "query_id",      # Query sequence ID
    "subject_id",    # Subject (database) sequence ID
    "percent_identity",  # Percentage of identical matches
    "alignment_length",  # Alignment length
    "mismatches",    # Number of mismatches
    "gap_opens",     # Number of gap openings
    "query_start",   # Start of alignment in query
    "query_end",     # End of alignment in query
    "subject_start", # Start of alignment in subject
    "subject_end",   # End of alignment in subject
    "evalue",        # Expectation value
    "bit_score"      # Bit score
    ]

    def parse_blast_result(result_text: str) -> pd.DataFrame:
        """
        Parses BLAST outfmt=6 result text into a pandas DataFrame.

        Args:
            result_text (str): BLAST result in outfmt=6.

        Returns:
            pd.DataFrame: BLAST results as a DataFrame.
        """
        # Split the result into lines and remove empty lines
        lines = [line for line in result_text.split("\n") if line.strip()]
        
        # Convert the lines into a DataFrame
        df = pd.DataFrame([line.split("\t") for line in lines], columns=BLAST_COLUMNS)
        
        # Convert numerical columns to appropriate types
        for col in ["percent_identity", "alignment_length", "mismatches", "gap_opens",
                    "query_start", "query_end", "subject_start", "subject_end",
                    "evalue", "bit_score"]:
            df[col] = pd.to_numeric(df[col])
        
        return df

    def blastp(self, query, db, evalue=0.001, outfmt=10, num_threads=50, max_target_seqs=50):
        """
        Perform a BLASTP search using the provided parameters.

        Args:
            query (str): The query sequence or filename.
            db (str): The database to search against.
            evalue (float): E-value threshold for BLAST hits. Default is 0.001.
            outfmt (int): Output format for BLAST results. Default is 10 (CSV).
            num_threads (int): Number of threads to use for the BLAST search. Default is 50.
            max_target_seqs (int): Maximum number of target sequences to retrieve. Default is 50.

        Returns:
            pd.DataFrame: A pandas DataFrame containing the BLAST results.
        """
        # now we will run the run service
        result = self.run_service(query, db, evalue, outfmt, num_threads, max_target_seqs)
        
        #parse the result to a pandas dataframe
        
        blast_result_text = result.text  # Assuming result.text contains the outfmt=6 text
        df = self.parse_blast_result(blast_result_text)
        
        return df

        
    
        
        