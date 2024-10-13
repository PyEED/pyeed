import asyncio
import io
import logging
from concurrent.futures import ThreadPoolExecutor
from enum import Enum, EnumMeta
from typing import List, Optional

from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Record import Blast as BlastRecord
from pydantic import BaseModel, Field

LOGGER = logging.getLogger(__name__)


class MetaEnum(EnumMeta):
    def __contains__(cls, item):
        try:
            cls(item)
        except ValueError:
            return False
        return True


class BaseEnum(Enum, metaclass=MetaEnum):
    pass


class BlastProgram(BaseEnum):
    BLASTP = "blastp"
    BLASTN = "blastn"
    BLASTX = "blastx"
    TBLASTN = "tblastn"
    TBLASTX = "tblastx"


class NCBIDataBase(BaseEnum):
    NR = "nr"
    UNIPROTKB = "swissprot"
    PDB = "pdb"
    REFSEQ = "refseq_protein"
    NT = "nt"


class SubstitutionMatrix(BaseEnum):
    BLOSUM45 = "BLOSUM45"
    BLOSUM62 = "BLOSUM62"
    BLOSUM80 = "BLOSUM80"
    PAM30 = "PAM30"
    PAM70 = "PAM70"


class Blast(BaseModel):
    query: str = Field(
        description="The query sequence",
        default=None,
    )

    n_hits: int = Field(
        description="Maximum number of hits to return",
        default=100,
    )

    evalue: float = Field(
        description="Expectation value (E) to safe hits",
        default=10,
    )

    matrix: str = Field(
        description="Substitution matrix",
        default=SubstitutionMatrix.BLOSUM62.value,
    )

    identity: float = Field(
        description="Minimum identity to accept hit",
        default=0.0,
        ge=0.0,
        le=1,
    )

    def run(self, program: str, ncbi_db: str) -> io.StringIO:
        assert (
            program in BlastProgram
        ), f"Invalid program: {program}, valid programs: {BlastProgram}"
        assert (
            ncbi_db in NCBIDataBase
        ), f"Invalid database: {ncbi_db}, valid databases: {NCBIDataBase}"

        if program == BlastProgram.BLASTP.value:
            return NCBIWWW.qblast(
                program,
                ncbi_db,
                self.query,
                expect=self.evalue,
                matrix_name=self.matrix,
                hitlist_size=self.n_hits,
            )

        elif program == BlastProgram.BLASTN.value:
            return NCBIWWW.qblast(
                program=program,
                database=ncbi_db,
                sequence=self.query,
                expect=self.evalue,
                hitlist_size=self.n_hits,
            )

    async def async_run(
        self,
        ncbi_db: str,
        program: str,
        foreign_executor: Optional[ThreadPoolExecutor] = None,
    ) -> io.StringIO:
        if not foreign_executor:
            executor = ThreadPoolExecutor()
        else:
            executor = foreign_executor

        loop = asyncio.get_running_loop()

        with executor as pool:
            result = await loop.run_in_executor(pool, self.run, program, ncbi_db)

        if not foreign_executor:
            executor.shutdown()

        return result

    def read(self, result: io.StringIO) -> BlastRecord:
        return NCBIXML.read(result)

    def extract_accession(self, record: io.StringIO) -> List[str]:
        record = NCBIXML.read(record)

        hits = []
        for hit in record.alignments:
            if hit.hsps[0].identities / hit.hsps[0].align_length > self.identity:
                hits.append(hit.accession)

        return hits


if __name__ == "__main__":
    from rich.status import Status

    async def main():
        seq = "MRNINVQLNPLSDIEKLQVELVERKGLGHPDYIADAVAEEASRKLSLYYLKKYGVILHHNLDKTLVVGGQATPRFKGGDVIQPIYIVVAGRATTEVKTESGIEQIPVGTIIIESVKEWIRNNFRYLDAEKHLIVDYKIGKGSTDLVGIFEAGKRVPLSNDTSFGVGFAPFTKLEKLVYETERHLNSKQFKAKLPEVGEDIKVMGLRRGNEVDLTIAMATISELIEDVNHYINVKEQAKNKILDLASKIAPDYDVRIYVNTGDKIDKNILYLTVTGTSAEHGDDGMTGRGNRGVGLITPMRPMSLEATAGKNPVNHVGKLYNVLANLIANKIAQEVKDVKFSQVQVLGQIGRPIDDPLIANVDVITYDGKLNDETKNEISGIVDEMLSSFNKLTELILEGKATLF"
        blast = Blast(query=seq, n_hits=10)
        executor = ThreadPoolExecutor(max_workers=4)

        blast_task = asyncio.create_task(
            blast.async_run(
                NCBIDataBase.UNIPROTKB.value, BlastProgram.BLASTP.value, executor
            )
        )

        with Status("Running BLAST"):
            result = await blast_task

        executor.shutdown()

        return result

    results = asyncio.run(main())

    print(results)
