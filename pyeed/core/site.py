from typing import Dict, List, Optional
from uuid import uuid4

from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from .abstractannotation import AbstractAnnotation


class Site(
    AbstractAnnotation,
    search_mode="unordered",
):
    """Position(s) constituting a site within a sequence."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    positions: List[int] = element(
        description="Position of the site(s) within the sequence.",
        default_factory=ListPlus,
        tag="positions",
        json_schema_extra=dict(
            multiple=True,
            term="http://edamontology.org/data_1016",
        ),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="d74df50f2e1144d97bc505140af6f320acba7e80"
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
