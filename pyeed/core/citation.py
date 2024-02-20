import sdRDM

from typing import Dict, List, Optional
from pydantic import PrivateAttr, model_validator
from uuid import uuid4
from pydantic_xml import attr, element, wrapped
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from .author import Author


@forge_signature
class Citation(
    sdRDM.DataModel,
    nsmap={
        "": "https://github.com/PyEED/pyeed@3b002efa6bf51e951767d8a7749ebad563897cb8#Citation"
    },
):
    """Information on publication of the entry 📖"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    doi: Optional[str] = element(
        description="DOI of the publication",
        default=None,
        tag="doi",
        json_schema_extra=dict(),
    )

    pubmed_id: Optional[str] = element(
        description="PubMed ID of the publication",
        default=None,
        tag="pubmed_id",
        json_schema_extra=dict(),
    )

    medline_id: Optional[str] = element(
        description="Medline ID of the publication",
        default=None,
        tag="medline_id",
        json_schema_extra=dict(),
    )

    year: Optional[int] = element(
        description="Year of publication",
        default=None,
        tag="year",
        json_schema_extra=dict(),
    )

    authors: List[Author] = wrapped(
        "authors",
        element(
            description="Authors of the publication",
            default_factory=ListPlus,
            tag="Author",
            json_schema_extra=dict(multiple=True),
        ),
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

    def add_to_authors(
        self,
        given_name: Optional[str] = None,
        family_name: Optional[str] = None,
        id: Optional[str] = None,
    ) -> Author:
        """
        This method adds an object of type 'Author' to attribute authors

        Args:
            id (str): Unique identifier of the 'Author' object. Defaults to 'None'.
            given_name (): Given name of the author. Defaults to None
            family_name (): Family name of the author. Defaults to None
        """
        params = {"given_name": given_name, "family_name": family_name}
        if id is not None:
            params["id"] = id
        self.authors.append(Author(**params))
        return self.authors[-1]
