import sdRDM

from typing import Dict, Optional
from pydantic import PrivateAttr, model_validator
from uuid import uuid4
from pydantic_xml import attr, element
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from .citation import Citation
from .organism import Organism


@forge_signature
class AbstractSequence(
    sdRDM.DataModel,
    nsmap={
        "": "https://github.com/PyEED/pyeed@3b002efa6bf51e951767d8a7749ebad563897cb8#AbstractSequence"
    },
):
    """"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    source_id: Optional[str] = element(
        description="Identifier of the sequence in the source database",
        default=None,
        tag="source_id",
        json_schema_extra=dict(),
    )

    name: Optional[str] = element(
        description="Name of the sequence",
        default=None,
        tag="name",
        json_schema_extra=dict(),
    )

    sequence: str = element(
        description="Sequence of the molecule",
        tag="sequence",
        json_schema_extra=dict(),
    )

    organism: Optional[Organism] = element(
        description="Corresponding organism",
        default=None,
        tag="organism",
        json_schema_extra=dict(),
    )

    citation: Optional[Citation] = element(
        description="Publication of the sequence",
        default_factory=Citation,
        tag="citation",
        json_schema_extra=dict(),
    )
    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="3b002efa6bf51e951767d8a7749ebad563897cb8"
    )
    _raw_xml_data: Dict = PrivateAttr(default_factory=dict)

    @model_validator(mode="after")
    def _parse_raw_xml_data(self):
        for attr, value in self:
            if isinstance(value, (ListPlus, list)) and all(
                (isinstance(i, _Element) for i in value)
            ):
                self._raw_xml_data[attr] = [elem2dict(i) for i in value]
            elif isinstance(value, _Element):
                self._raw_xml_data[attr] = elem2dict(value)
        return self

    def _fasta_string(self):
        return f">{self.source_id}\n{self.sequence}"

    def to_fasta(self):
        return self._fasta_string()
