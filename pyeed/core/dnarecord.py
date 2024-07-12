from typing import Dict, Optional, List
from uuid import uuid4

from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from IPython.display import clear_output
from sdRDM.tools.utils import elem2dict
import asyncio

from pyeed.fetch.dnafetcher import DNAFetcher


from pyeed.core.sequencerecord import SequenceRecord


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

        return asyncio.run(
            DNAFetcher(ids=accession_ids).fetch(force_terminal=False)
        )
    

if __name__ == "__main__":
    dna_record = DNARecord.get_id('AF188200.1')

    print(f"DNA Record: {dna_record}")

    