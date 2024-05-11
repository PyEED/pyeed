import io
from typing import Dict, List

import httpx
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from pyeed.tools.abstract_tool import AbstractTool, ServiceURL


class ClustalOmega(AbstractTool):
    """
    Class for ClustalOmega aligner running as a REST service.
    """

    _service_url = ServiceURL.CLUSTALO.value

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
        file = self.create_file(data)
        return httpx.post(self._service_url, files=file, timeout=600)

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
