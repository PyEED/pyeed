from typing import List, Dict, Any, Optional

import httpx
from pydantic import PrivateAttr, BaseModel, Field

from pyeed.tools.abstract_tool import AbstractTool, ServiceURL

## IN AND OUT DATA MODELS -------------

class BlastRequest(BaseModel):
    tool: Optional[str] = 'blastp'
    query: str
    db: str
    evalue: Optional[str] = '0.001'
    outfmt: Optional[str] = '6'

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
        
    def run_service(self):
        """Executes the service."""
        pass


if __name__ == "__main__":

    blastp = BlastP()

    # Check if the service is running
    response = blastp.check_running_service()
    print(response.status_code)