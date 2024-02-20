import sdRDM

from typing import Dict, Optional
from pydantic import PrivateAttr, model_validator
from uuid import uuid4
from pydantic_xml import attr, element
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from .citation import Citation


@forge_signature
class Substrate(
    sdRDM.DataModel,
    nsmap={
        "": "https://github.com/PyEED/pyeed@3b002efa6bf51e951767d8a7749ebad563897cb8#Substrate"
    },
):
    """Promiscuous substrate of an enzyme 🧪"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    name: Optional[str] = element(
        description="Name of the substrate",
        default=None,
        tag="name",
        json_schema_extra=dict(),
    )

    inchi: Optional[str] = element(
        description="InChI code of the substrate",
        default=None,
        tag="inchi",
        json_schema_extra=dict(),
    )

    smiles: Optional[str] = element(
        description="SMILES code of the substrate",
        default=None,
        tag="smiles",
        json_schema_extra=dict(),
    )

    chebi_id: Optional[str] = element(
        description="ChEBI ID of the substrate",
        default=None,
        tag="chebi_id",
        json_schema_extra=dict(),
    )

    citation: Optional[Citation] = element(
        description="Citations of the substrate",
        default_factory=Citation,
        tag="citation",
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
