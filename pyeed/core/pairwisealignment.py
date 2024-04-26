from typing import Dict, Optional
from uuid import uuid4

from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from .alignmentdata import AlignmentData


class PairwiseAlignment(
    AlignmentData,
    search_mode="unordered",
):
    """"""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    score: Optional[float] = element(
        description="Alignment score",
        default=None,
        tag="score",
        json_schema_extra=dict(),
    )

    identity: Optional[float] = element(
        description="Ration of identical residues in the alignment",
        default=None,
        tag="identity",
        json_schema_extra=dict(),
    )

    similarity: Optional[float] = element(
        description="Ration of similar residues in the alignment",
        default=None,
        tag="similarity",
        json_schema_extra=dict(),
    )

    gaps: Optional[int] = element(
        description="Number of gaps in the alignment",
        default=None,
        tag="gaps",
        json_schema_extra=dict(),
    )

    mismatches: Optional[int] = element(
        description="Number of mismatches in the alignment",
        default=None,
        tag="mismatches",
        json_schema_extra=dict(),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="935e9b5cb924657615539e6d4813287fea72ff35"
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
