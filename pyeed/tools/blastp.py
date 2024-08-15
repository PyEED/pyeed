import io
import json
from typing import List, Dict, Any, Optional

import httpx
from Bio.SearchIO import BlastIO
from pydantic import PrivateAttr, BaseModel, Field

from pyeed.tools.abstract_tool import AbstractTool, ServiceURL

class BlastP(AbstractTool):
    """
    Class for BLASTP running as a REST service.
    """

    _service_url = PrivateAttr(ServiceURL.BLAST_LOCAL.value)

    def check_running_service(self):
        # Check if the service is running, just check the status code
        try:
            return httpx.put(self._service_url)

        except httpx.ConnectError as connect_error:
            error_number = connect_error.__context__.args[0].errno
            if error_number == 8 or error_number == -3:
                self._service_url = self._service_url.replace("blast_docker", "localhost")
                try:
                    return httpx.get(self._service_url)
                except httpx.ConnectError as connect_error:
                    raise httpx.ConnectError("Blast Docker Service is not running.")

            raise httpx.ConnectError("Blast Docker Service is not running.")
        
    def extact_ids_from_blast_string(self, response):
        # here we extract the result from the written output file
        # the file is out.out as specified in FASTAPI for our BLASTP service
        data = response.json()
        import Bio.SearchIO as SearchIO
        
        blast_records = SearchIO.read(io.StringIO(data), 'blast-tab')

        ids = [ hit.id for hit in blast_records ]

        return ids

    
        
    def blastp(self, query: str, db: str, evalue: str, outfmt: str, num_threads: str):
        # here we run the actual search
        data = {
            "query": query,
            "db": db,
            "evalue": evalue,
            "outfmt": outfmt,
            "tool": "blastp",
            "num_threads": num_threads
        }

        result = self.run_service(data)

        if result.status_code == 500:
            raise httpx.HTTPStatusError("Error in running BLASTP service. The error message is: " + result.text)
        # now we need to extract the result from the written output file
        # before that we need to ensure that the service is already finished
        return self.extact_ids_from_blast_string(result)

        

    def run_service(self, data) -> httpx.Response:
        # here we run the actual search
        try:
            return httpx.post(self._service_url, timeout=6000, json=data, headers = {'Content-Type': 'application/json'})

        except httpx.ConnectError as connect_error:
            error_number = connect_error.__context__.args[0].errno
            if error_number == 8 or error_number == -3:
                self._service_url = self._service_url.replace("blast_docker", "localhost")
                try:
                    return httpx.post(self._service_url, json=data, timeout=6000, headers = {'Content-Type': 'application/json'})
                except httpx.ConnectError as connect_error:
                    raise httpx.ConnectError("Blast Docker Service is not running.")

            raise httpx.ConnectError("Blast Docker Service is not running.")



if __name__ == "__main__":

    from pyeed.core import ProteinRecord

    mat_accessions = [
    "MBP1912539.1",
    "SEV92896.1",
    "MBO8174569.1",
    "WP_042680787.1",
    "NPA47376.1",
    "WP_167889085.1",
    "WP_048165429.1",
    "ACS90033.1",
    ]

    mats = ProteinRecord.get_ids(mat_accessions)

    ids = mats[0].ncbi_blast_local(db = '/blast/blastdb/nr/nr', evalue=0.005, num_threads=50)

    print(ids)