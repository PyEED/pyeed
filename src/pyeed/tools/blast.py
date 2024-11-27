import io

import httpx
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Record import Blast as BlastRecord
from pyeed.tools.abstract_tool import AbstractTool, ServiceURL

class Blast(AbstractTool):
    """
        Class for Blast search with biopython
    """

def __init__(self):
    super().__init__()
    self._service_url = ServiceURL.BLAST.value

def create_file(self, multifasta: list[str]) -> dict[str, str]:
    """
        Sets up the input data for the Blast container.

        Args:
            multifasta (list[str]): List of FASTA formatted sequences to be aligned.

        Returns:
            dict[str, str]: A dictionary containing the input file.
    """
    data = "\n".join(multifasta)
    return {"file": data}

def run_service(self, data) --> httpx.Response:
    """
        Executes the Blast service.
    """
    file = self.create_file(data)
    try: 
        return httpx.post(self._service_url, files=file, timeout=600)
    
    except httpx.ConnectError as connect_error:
        context = connect_error.__context__
        if context and hasattr(context, "args"):
            error_number = context.args[0].errno
            if error_number == 8 or error_number == -3:
                self._service_url = self._service_url.replace(
                    "blast", "localhost"
                )
                try:
                    return httpx.post(self._service_url, files=file, timeout=600)
                except httpx.ConnectError:
                    raise httpx.ConnectError(
                        "PyEED Docker Service is not running."
                    ) from None
        
        print(connect_error)
        raise httpx.ConnectError("PyEED Docker Service is not running.") from None
    
    
def blast(
    self, 
    sequence: str, 
    database: str = "nr",
    substitution_matrix: str = "BLOSSUM62", 
    program: str = "blastp"
    ) -> BlastRecord: 
    """
        Runs a Blast search with the given sequence.

        Args:
            sequence (str): The sequence to search for.

        Returns:
            BlastRecord: The Blast search result.
    """
    
    BLASTP = "blastp"
    BLASTN = "blastn"
    BLASTX = "blastx"
    TBLASTN = "tblastn"
    TBLASTX = "tblastx"


    NR = "nr"
    UNIPROTKB = "swissprot"
    PDB = "pdb"
    REFSEQ = "refseq_protein"
    NT = "nt"

    BLOSUM45 = "BLOSUM45"
    BLOSUM62 = "BLOSUM62"
    BLOSUM80 = "BLOSUM80"
    PAM30 = "PAM30"
    PAM70 = "PAM70"
    
    if program == "blastp":
        return NCBIWWW.qblast(
            program=program,
            database=database,
            sequence=sequence,
            expect=self.evalue,
            hitlist_size=self.n_hits,
        )
    
    result = NCBIWWW.qblast("blastp", "nr", sequence)
    return NCBIXML.read(io.StringIO(result))