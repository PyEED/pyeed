import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .author import Author


@forge_signature
class Citation(sdRDM.DataModel):
    """Information on publication of the entry 📖"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("citationINDEX"),
        xml="@id",
    )

    doi: Optional[str] = Field(
        default=None,
        description="DOI of the publication",
    )

    pubmed_id: Optional[str] = Field(
        default=None,
        description="PubMed ID of the publication",
    )

    medline_id: Optional[str] = Field(
        default=None,
        description="Medline ID of the publication",
    )

    year: Optional[int] = Field(
        default=None,
        description="Year of publication",
    )

    authors: List[Author] = Field(
        description="Authors of the publication",
        default_factory=ListPlus,
        multiple=True,
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    __commit__: Optional[str] = PrivateAttr(
        default="2c478e9b9618bfdc095c0c8906fbe67c80a3e2d7"
    )

    def add_to_authors(
        self,
        given_name: Optional[str] = None,
        family_name: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
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
