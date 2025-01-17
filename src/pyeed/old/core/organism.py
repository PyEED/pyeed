from typing import Dict, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict


class Organism(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """Description of an organism 🦠."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    taxonomy_id: int = element(
        description=(
            "A stable unique identifier for each taxon (for a species, a family, an"
            " order, or any other group in the NCBI taxonomy database."
        ),
        tag="taxonomy_id",
        json_schema_extra=dict(
            term="http://edamontology.org/data_1179",
        ),
    )

    name: Optional[str] = element(
        description="The name of an organism (or group of organisms).",
        default=None,
        tag="name",
        json_schema_extra=dict(
            term="http://edamontology.org/data_2909",
        ),
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
        json_schema_extra=dict(
            term="http://edamontology.org/data_1044",
        ),
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
        description="The name of a family of organism.",
        default=None,
        tag="family",
        json_schema_extra=dict(
            term="http://edamontology.org/data_2732",
        ),
    )

    genus: Optional[str] = element(
        description="The name of a genus of organism.",
        default=None,
        tag="genus",
        json_schema_extra=dict(
            term="http://edamontology.org/data_1870",
        ),
    )

    species: Optional[str] = element(
        description="The name of a species (typically a taxonomic group) of organism.",
        default=None,
        tag="species",
        json_schema_extra=dict(
            term="http://edamontology.org/data_1045",
        ),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="72d2203f2e3ce4b319b29fa0d2f146b5eead7b00"
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
