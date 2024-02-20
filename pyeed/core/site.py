import sdRDM

from typing import Dict, List, Optional
from pydantic import PrivateAttr, model_validator
from uuid import uuid4
from pydantic_xml import attr, element, wrapped
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from .proteinsitetype import ProteinSiteType


@forge_signature
class Site(
    sdRDM.DataModel,
    nsmap={
        "": "https://github.com/PyEED/pyeed@3b002efa6bf51e951767d8a7749ebad563897cb8#Site"
    },
):
    """Annotation of a site within a sequence 📍"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    name: Optional[str] = element(
        description="Name of the site",
        default=None,
        tag="name",
        json_schema_extra=dict(),
    )

    type: Optional[ProteinSiteType] = element(
        description="Type of the site",
        default=None,
        tag="type",
        json_schema_extra=dict(),
    )

    positions: List[int] = wrapped(
        "positions",
        element(
            description="Positions of the site",
            default_factory=ListPlus,
            tag="integer",
            json_schema_extra=dict(multiple=True),
        ),
    )

    cross_ref: Optional[str] = element(
        description="Database cross reference",
        default=None,
        tag="cross_ref",
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
