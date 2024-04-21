from typing import Dict, List, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict


class StandardNumbering(
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

    reference_accession_id: Optional[str] = element(
        description="Standard numbering of the reference sequence",
        default=None,
        tag="reference_accession_id",
        json_schema_extra=dict(),
    )

    numbered_accession_id: Optional[str] = element(
        description="Standard numbering of the query sequence",
        default=None,
        tag="numbered_accession_id",
        json_schema_extra=dict(),
    )

    numbering: List[str] = element(
        description="Standard numbering of the aligned sequence",
        default_factory=ListPlus,
        tag="numbering",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="9cf5a53e2161d63c639e27dadd2d0bda5c174f67"
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
