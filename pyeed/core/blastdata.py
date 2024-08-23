from typing import Dict, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict


class BlastData(
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

    identity: Optional[float] = element(
        description="Minimum identity to safe hits.",
        default=0.0,
        tag="identity",
        json_schema_extra=dict(),
    )

    evalue: Optional[float] = element(
        description="Expectation value (E) to safe hits.",
        default=10.0,
        tag="evalue",
        json_schema_extra=dict(),
    )

    n_hits: Optional[int] = element(
        description="Number of hits to return.",
        default=100,
        tag="n_hits",
        json_schema_extra=dict(),
    )

    substitution_matrix: Optional[str] = element(
        description="Substitution matrix to use.",
        default="'BLOSUM62'",
        tag="substitution_matrix",
        json_schema_extra=dict(),
    )

    word_size: Optional[int] = element(
        description="Word size of the initial match.",
        default=3,
        tag="word_size",
        json_schema_extra=dict(
            inclusivminimum=2,
            inclusivemaximum=7,
        ),
    )

    gap_open: Optional[float] = element(
        description="Gap open penalty.",
        default=11.0,
        tag="gap_open",
        json_schema_extra=dict(),
    )

    gap_extend: Optional[float] = element(
        description="Gap extend penalty.",
        default=1.0,
        tag="gap_extend",
        json_schema_extra=dict(),
    )

    threshold: Optional[float] = element(
        description="Minimum score to add a word to the BLAST lookup table.",
        default=11,
        tag="threshold",
        json_schema_extra=dict(),
    )

    db_name: Optional[str] = element(
        description="Name of the database to search.",
        default=None,
        tag="db_name",
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
