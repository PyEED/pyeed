from typing import Dict, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from .organism import Organism


class SequenceRecord(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """A molecular sequence and associated annotation data."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    uri: Optional[str] = element(
        description="URI of the sequence.",
        default=None,
        tag="uri",
        json_schema_extra=dict(
            term="http://edamontology.org/data_1047",
        ),
    )

    accession_id: Optional[str] = element(
        description="Accession ID of the sequence.",
        default=None,
        tag="accession_id",
        json_schema_extra=dict(
            term="http://edamontology.org/data_2091",
        ),
    )

    name: Optional[str] = element(
        description="Arbitrary name of the sequence.",
        default=None,
        tag="name",
        json_schema_extra=dict(
            term="http://edamontology.org/data_2099",
        ),
    )

    organism: Optional[Organism] = element(
        description="The organism from which the sequence was obtained.",
        default=None,
        tag="organism",
        json_schema_extra=dict(
            term="http://edamontology.org/data_2530",
        ),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="8fa955f44bf4de5d3aff5efeed401995ab415498"
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
