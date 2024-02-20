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
class Organism(
    sdRDM.DataModel,
    nsmap={
        "": "https://github.com/PyEED/pyeed@3b002efa6bf51e951767d8a7749ebad563897cb8#Organism"
    },
):
    """Description of an organism 🦠"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    name: Optional[str] = element(
        description="Name of the organism",
        default=None,
        tag="name",
        json_schema_extra=dict(),
    )

    taxonomy_id: str = element(
        description="NCBI Taxonomy ID to identify the organism",
        tag="taxonomy_id",
        json_schema_extra=dict(),
    )

    domain: Optional[str] = element(
        description="Domain of the organism",
        default=None,
        tag="domain",
        json_schema_extra=dict(),
    )

    kingdom: Optional[str] = element(
        description="Kingdom of the organism",
        default=None,
        tag="kingdom",
        json_schema_extra=dict(),
    )

    phylum: Optional[str] = element(
        description="Phylum of the organism",
        default=None,
        tag="phylum",
        json_schema_extra=dict(),
    )

    tax_class: Optional[str] = element(
        description="Class of the organism",
        default=None,
        tag="tax_class",
        json_schema_extra=dict(),
    )

    order: Optional[str] = element(
        description="Order of the organism",
        default=None,
        tag="order",
        json_schema_extra=dict(),
    )

    family: Optional[str] = element(
        description="Family of the organism",
        default=None,
        tag="family",
        json_schema_extra=dict(),
    )

    genus: Optional[str] = element(
        description="Genus of the organism",
        default=None,
        tag="genus",
        json_schema_extra=dict(),
    )

    species: Optional[str] = element(
        description="Species of the organism",
        default=None,
        tag="species",
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
