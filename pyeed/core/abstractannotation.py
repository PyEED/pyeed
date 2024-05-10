from typing import Dict, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict


class AbstractAnnotation(
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

    uri: Optional[str] = element(
        description="URI of the annotation.",
        default=None,
        tag="uri",
        json_schema_extra=dict(
            term="http://edamontology.org/data_1047",
        ),
    )

    accession_id: Optional[str] = element(
        description="Accession ID of the annotation.",
        default=None,
        tag="accession_id",
        json_schema_extra=dict(
            term="http://edamontology.org/data_2091",
        ),
    )

    name: Optional[str] = element(
        description="A name of a sequence feature, e.g. the name of a feature",
        default=None,
        tag="name",
        json_schema_extra=dict(
            term="http://edamontology.org/data_2099",
        ),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="a7defc5c87a2296a2e4b522b07236b2aef6413ac"
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
