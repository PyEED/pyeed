import io
import json
from typing import List, Dict, Any, Optional

import httpx
from Bio.SearchIO import BlastIO
from pydantic import PrivateAttr, BaseModel, Field

from pyeed.tools.abstract_tool import AbstractTool, ServiceURL

class MMSeqs(AbstractTool):

    _service_url = PrivateAttr(ServiceURL.MMSEQS_LOCAL.value)

    def check_running_service(self):
        """
        Check if the service is running, just check the status code

        @app.get("/")
        async def read_root():
            return {"message": "Welcome to the MMSeqs2 API!"}

        Returns:
            [type]: [description]
        """
        try:
            return httpx.get(self._service_url).text == "Welcome to the MMSeqs2 API!"
        
        except httpx.ConnectError as connect_error:
            error_number = connect_error.__context__.args[0].errno
            if error_number == 8 or error_number == -3:
                self._service_url = self._service_url.replace("mmseqs_docker", "localhost")
                try:
                    return httpx.get(self._service_url).text == "Welcome to the MMSeqs2 API!"
                except httpx.ConnectError as connect_error:
                    raise httpx.ConnectError("MMSeqs Docker Service is not running.")

            raise httpx.ConnectError("MMSeqs Docker Service is not running.")
        
    def extact_ids_from_results_file(self, job_id):
        try:
            response = httpx.get(f"{self._service_url}results/{job_id}")

        except httpx.ConnectError as connect_error:
            error_number = connect_error.__context__.args[0].errno
            if error_number == 8 or error_number == -3:
                self._service_url = self._service_url.replace("mmseqs_docker", "localhost")
                try:
                    response = httpx.get(f"{self._service_url}results/{job_id}")
                except httpx.ConnectError as connect_error:
                    raise httpx.ConnectError("MMSeqs Docker Service is not running.")

            raise httpx.ConnectError("MMSeqs Docker Service is not running.")
        
        data = response.json()

        import Bio.SearchIO as SearchIO
        
        blast_records = SearchIO.read(io.StringIO(data), 'blast-tab')

        ids = [ hit.id for hit in blast_records ]

        return ids
            

    def run_service(self, data):
        try:
            return httpx.post(self._service_url, timeout=6000, json=data, headers = {'Content-Type': 'application/json'})

        except httpx.ConnectError as connect_error:
            error_number = connect_error.__context__.args[0].errno
            if error_number == 8 or error_number == -3:
                self._service_url = self._service_url.replace("mmseqs_docker", "localhost")
                try:
                    return httpx.post(self._service_url, json=data, timeout=6000, headers = {'Content-Type': 'application/json'})
                except httpx.ConnectError as connect_error:
                    raise httpx.ConnectError("MMSeqs Docker Service is not running.")

            raise httpx.ConnectError("MMSeqs Docker Service is not running.")

    def mmseqs(self, query: str, database: str, num_threads: str, sensitivity: str):
        data = {
            "query": query,
            "database": database,
            "num_threads": num_threads,
            "sensitivity": sensitivity,
            "blast_format": True
            }

        job_id = self.run_service(data)

        results_biopython = self.extact_ids_from_results_file(job_id)


        return results_biopython
        



if __name__ == "__main__":

    from pyeed.core import ProteinRecord

    mmseqs = MMSeqs()
    
    check = mmseqs.check_running_service()
    print(check)

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

    mmseqs = MMSeqs()

    ids = mmseqs.mmseqs(mats[0].seq, "nr_mmseqs", "40", "5.7")

    print(ids)

