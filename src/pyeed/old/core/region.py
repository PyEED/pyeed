from typing import Dict, Optional
from uuid import uuid4

from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from .abstractannotation import AbstractAnnotation


class Region(
    AbstractAnnotation,
    search_mode="unordered",
):
    """Regional annotation of a feature within a sequence. The direction of the region is defined by the start and end positions."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    start: Optional[int] = element(
        description="Start position of the site.",
        default=None,
        tag="start",
        json_schema_extra=dict(
            term="http://semanticscience.org/resource/SIO_000943",
        ),
    )

    end: Optional[int] = element(
        description="End position of the site.",
        default=None,
        tag="end",
        json_schema_extra=dict(
            term="http://semanticscience.org/resource/SIO_000953",
        ),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="b926bfec3aa1ec45a5614cf6ac4a546252dd384c"
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
