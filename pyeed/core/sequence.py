from typing import Dict, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict


class Sequence(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """"""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    sequence_id: Optional[str] = element(
        description="Identifier of the sequence in the source database",
        default=None,
        tag="sequence_id",
        json_schema_extra=dict(),
    )

    sequence: Optional[str] = element(
        description="Molecular sequence.",
        default=None,
        tag="sequence",
        json_schema_extra=dict(),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="eced817a915618922cb780cdd0025d52b04b159d"
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

    def fasta_string(self) -> str:
        """
        This method returns the sequence in FASTA format

        Returns:
            str: Sequence in FASTA format
        """
        return f">{self.source_id}\n{self.sequence}"

    def __str__(self) -> str:
        return self.fasta_string()
