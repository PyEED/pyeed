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
        
        #save the result as a csv file and read it with pandas

        
        with open('blast_result.csv', 'w') as f:
            f.write(result.json()['result'])
            
        # reed the csv file
        df = pd.read_csv('blast_result.csv')
        
        return df.head()

        
    
        
        