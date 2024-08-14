import json
from typing import List, Dict, Any, Optional

import httpx
from pydantic import PrivateAttr, BaseModel, Field

from pyeed.tools.abstract_tool import AbstractTool, ServiceURL

class BlastP(AbstractTool):
    """
    Class for BLASTN running as a REST service.
    """

    _service_url = PrivateAttr(ServiceURL.BLAST_LOCAL.value)

    def check_running_service(self):
        # Check if the service is running
        try:
            return httpx.get(self._service_url)

        except httpx.ConnectError as connect_error:
            error_number = connect_error.__context__.args[0].errno
            if error_number == 8 or error_number == -3:
                self._service_url = self._service_url.replace("blast_docker", "localhost")
                try:
                    return httpx.get(self._service_url)
                except httpx.ConnectError as connect_error:
                    raise httpx.ConnectError("Blast Docker Service is not running.")

            raise httpx.ConnectError("Blast Docker Service is not running.")
        
    def blastp(self, query: str, db: str, evalue: str, outfmt: str):
        # here we run the actual search
        data = {
            "query": query,
            "db": db,
            "evalue": evalue,
            "outfmt": outfmt,
            "tool": "blastp"
        }

        result = self.run_service(data)
        print(result)
        print(result.read())
        

    def run_service(self, data) -> httpx.Response:
        # here we run the actual search
        try:
            print(data)
            return httpx.post(self._service_url, timeout=600, json=data, headers = {'Content-Type': 'application/json'})

        except httpx.ConnectError as connect_error:
            error_number = connect_error.__context__.args[0].errno
            if error_number == 8 or error_number == -3:
                self._service_url = self._service_url.replace("blast_docker", "localhost")
                try:
                    return httpx.post(self._service_url, json=data, timeout=600, headers = {'Content-Type': 'application/json'})
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

    blastp = BlastP()

    # Check if the service is running
    response = blastp.check_running_service()
    print(response.status_code)

    print(mats[0].sequence)

    blastp.blastp(query=mats[0].sequence, db="nr", evalue="0.001", outfmt="6")