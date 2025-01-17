import io
from typing import Dict, List

import httpx
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from pydantic import PrivateAttr
from pyeed.tools.abstract_tool import AbstractTool, ServiceURL


class ClustalOmega(AbstractTool):
    """
    Class for ClustalOmega aligner running as a REST service.
    """

    _service_url = PrivateAttr(ServiceURL.CLUSTALO.value)

    def create_file(self, multifasta: List[str]) -> Dict[str, str]:
        """
        Sets up the input data for the ClustalOmega container.

        Args:
            multifasta (List[str]): List of FASTA formatted sequences to be aligned.

        Returns:
            Dict[str, str]: A dictionary containing the input file.
        """
        data = "\n".join(multifasta)
        return {"file": data}

    def extract_output_data(self, response) -> MultipleSeqAlignment:
        """
        Extracts the output data from the ClustalOmega container.

        Returns:
            MultiSequenceAlignment: The alignment result.
        """
        alignment = AlignIO.read(io.StringIO(response), "clustal")
        self._delete_temp_dir()

        return alignment

    def run_service(self, data) -> httpx.Response:
        """Executes the ClustalOmega service."""
        file = self.create_file(data)
        try:
            return httpx.post(self._service_url, files=file, timeout=600)

        except httpx.ConnectError as connect_error:
            error_number = connect_error.__context__.args[0].errno
            if error_number == 8 or error_number == -3:
                self._service_url = self._service_url.replace("clustalo", "localhost")
                try:
                    return httpx.post(self._service_url, files=file, timeout=600)
                except httpx.ConnectError as connect_error:
                    raise httpx.ConnectError("PyEED Docker Service is not running.")
            print(connect_error)
            raise httpx.ConnectError("PyEED Docker Service is not running.")

    def align(self, sequences: List[str]):
        """
        Aligns multiple sequences and returns the alignment result.

        Args:
            sequences (List[str]): List of FASTA formatted sequences to be aligned.

        Returns:
            MultiSequenceAlignment: The alignment result.
        """
        r = self.run_service(sequences)
        cleaned_text = r.text.replace('"', "").encode().decode("unicode_escape")

        return self.extract_output_data(cleaned_text)


if __name__ == "__main__":
    sequences = [
        ">seq1\nMTHKLLLTLLFTLLFSSAYSRG",
        ">seq2\nMTHKILLLTLLFTLLFSSAYSRG",
        ">seq3\nMTHKILLLTLLFTLLFSSCYSRG",
    ]

    clustalo = ClustalOmega()
    alignment = clustalo.align(sequences)
    print(alignment)
