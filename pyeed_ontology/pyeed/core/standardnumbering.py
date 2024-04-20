from typing import Dict, List, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict


@forge_signature
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
        default="f67fe20dfaf0d9cfcf4455c132f9f6d0bce52c06"
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