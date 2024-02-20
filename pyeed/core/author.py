import sdRDM

from typing import Dict, Optional
from pydantic import PrivateAttr, model_validator
from uuid import uuid4
from pydantic_xml import attr, element
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict


@forge_signature
class Author(
    sdRDM.DataModel,
    nsmap={
        "": "https://github.com/PyEED/pyeed@3b002efa6bf51e951767d8a7749ebad563897cb8#Author"
    },
):
    """"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    given_name: Optional[str] = element(
        description="Given name of the author",
        default=None,
        tag="given_name",
        json_schema_extra=dict(),
    )

    family_name: Optional[str] = element(
        description="Family name of the author",
        default=None,
        tag="family_name",
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
