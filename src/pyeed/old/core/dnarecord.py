import asyncio
import warnings
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional
from uuid import uuid4

from IPython.display import clear_output
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from pyeed.core.sequencerecord import SequenceRecord
from pyeed.fetch.blast import Blast, BlastProgram
from pyeed.fetch.dnafetcher import DNAFetcher
from rich.console import Console
from rich.status import Status
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict


class DNARecord(
    SequenceRecord,
    search_mode="unordered",
):
    """A nucleic acid sequence and associated metadata 🧬"""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    gc_content: Optional[float] = element(
        description="GC content of the sequence.",
        default=None,
        tag="gc_content",
        json_schema_extra=dict(),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="72d2203f2e3ce4b319b29fa0d2f146b5eead7b00"
    )

    _raw_xml_data: Dict = PrivateAttr(default_factory=dict)

    @model_validator(mode="after")
    def _parse_raw_xml_data(self):
        for attr, value in self:
            if isinstance(value, (ListPlus, list)) and all(
                isinstance(i, _Element) for i in value
            ):
                self._raw_xml_data[attr] = [elem2dict(i) for i in value]
            elif isinstance(value, _Element):
                self._raw_xml_data[attr] = elem2dict(value)

        return self

    @classmethod
    def get_id(cls, dna_id: str) -> "DNARecord":
        """
        This method creates a 'DNARecord' object from a given DNA accession ID.

        Args:
            dna_id (str): ID of the DNA in NCBI database.

        Returns:
            DNARecord: 'DNARecord' with information of the corresponding dna_id.
        """

        import nest_asyncio

        nest_asyncio.apply()

        if isinstance(dna_id, list) and all(isinstance(x, str) for x in dna_id):
            warnings.warn("For getting multiple sequences by ID use `get_ids` instead.")
            return cls.get_ids(dna_id)

        sequences = asyncio.run(DNAFetcher(ids=[dna_id]).fetch(quiet=False))[0]
        clear_output()
        return sequences

    @classmethod
    def get_ids(cls, accession_ids: List[str]) -> List["DNARecord"]:
        """Creates a list of 'DNARecord' objects from a list of dna accession IDs.

        Returns:
            List[DNARecord]: A list of 'DNARecord' objects representing the dna sequences found.
        """

        import nest_asyncio

        nest_asyncio.apply()

        return asyncio.run(DNAFetcher(ids=accession_ids).fetch(force_terminal=False))

    def ncbi_blast(
        self,
        n_hits: int,
        e_value: float = 10.0,
        db: str = "nt",
        identity: float = 0.0,
        **kwargs,
    ) -> List["DNARecord"]:
        """
        Runs a BLAST search using the NCBI BLAST service to find similar dna sequences.

        Args:
            n_hits (int): The number of hits to retrieve.
            e_value (float, optional): The maximum E-value threshold for reporting hits. Defaults to 10.0.
            db (str, optional): The database to search against. Defaults to "nt".
            match/mismatch (int, optional): Match/mismatch score. Defaults to 1/-2.
            gap_cost (int, optional): Cost to open a gap. Default is 0 and means linear
            identity (float, optional): The minimum sequence identity threshold for reporting hits. Defaults to 0.0.
            **kwargs: Additional keyword arguments.

        Returns:
            List[DNARecord]: A list of DNARecord objects representing the similar dna sequences found.

        Raises:
            AssertionError: If the specified database is not supported.

        Example:
            dna_info = DNARecord()
            similar_proteins = dna_info.ncbi_blast(n_hits=10, e_value=0.001, database="nt")
        """

        import nest_asyncio
        from pyeed.fetch.blast import NCBIDataBase

        nest_asyncio.apply()

        assert (
            db in NCBIDataBase
        ), f"Database needs to be one of {NCBIDataBase.__members__.keys()}"

        program = BlastProgram.BLASTN.value
        executor = ThreadPoolExecutor(max_workers=1)
        blaster = Blast(
            query=self.sequence,
            n_hits=n_hits,
            evalue=e_value,
            identity=identity,
        )

        with Status(
            "Running BLAST", console=Console(force_terminal=False, force_jupyter=False)
        ):
            result = asyncio.run(blaster.async_run(db, program, executor))
            clear_output()

        # result = blaster.run(program, db)

        accessions = blaster.extract_accession(result)

        return self.get_ids(accessions)


if __name__ == "__main__":
    dna_record = DNARecord.get_id("AF188200.1")

    print(f"DNA Record: {dna_record}")

    print("Running blast...")
    similar_dna_records = dna_record.ncbi_blast(n_hits=100, db="nt")
    print(similar_dna_records)
