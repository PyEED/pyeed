import io
import json
from typing import List, Dict, Any, Optional

import httpx
from Bio.SearchIO import BlastIO
from pydantic import PrivateAttr, BaseModel, Field

from pyeed.tools.abstract_tool import AbstractTool, ServiceURL

class BlastN(AbstractTool):
    """
    Class for BLASTN running as a REST service.

    Args:
        AbstractTool ([type]): [description]
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
        
    def extact_ids_from_blast_string(self, response) -> Dict[str, Any]:
        # here we extract the result from the written output file
        # the file is out.out as specified in FASTAPI for our BLASTP service
        data = response.json()
        import Bio.SearchIO as SearchIO
        
        blast_records = SearchIO.read(io.StringIO(data), 'blast-tab')

        ids = [ hit.id for hit in blast_records ]

        return ids

    def blastn(self, query: str, db: str, evalue: str, outfmt: str):
        # here we run the actual search
        data = {
            "query": query,
            "db": db,
            "evalue": evalue,
            "outfmt": outfmt,
            "tool": "blastn"
        }

        result = self.run_service(data)
        # now we need to extract the result from the written output file
        # before that we need to ensure that the service is already finished
        return self.extact_ids_from_blast_string(result)


    def run_service(self, data) -> httpx.Response:
        # here we run the actual search
        try:
            print(data)
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

    from pyeed.core import DNARecord

    # create single DNARecord
    dna = DNARecord(sequence="""ATGAGTATTCAACAGGATTGTCGTCGACCTGCGTGCCGGTGACCTGCAAGGCGCAGT
CAGTGCCGCGTGTTGACGTCGCGGCGAGGATGACCTGACCTGGTTGACCGGCAACCGT
TCCAGCAGTGAGCGGCCGAGGTGCCGAGCCAGTCCCTGACCGGCTGTTGATGTCGCC
CTGCCGACCTGGACGGTGCGGACCCGTCAGCAGGACCAGCTGACGTCGGGAGCAGT
AAGACCGCTGTCGCCGACGAGGAGGACGAGGACCTGTCGAGGACCTGCCGAGGACC""", id="test")

    blastp = BlastN()

    ids = blastp.blastp(query=dna.sequence, db='/blast/blastdb/nt/nt', evalue="0.001", outfmt="6")
    print(ids)