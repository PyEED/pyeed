from typing import Dict, List, Optional, Set
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from .region import Region


class RegionSet(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """A set of regions forming a higher order structure. For example, a set of exons in a gene, or a set of secondary structures forming a super-secondary structure."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    regions: List[Region] = element(
        description="Regions of the cluster.",
        default_factory=ListPlus,
        tag="regions",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="b926bfec3aa1ec45a5614cf6ac4a546252dd384c"
    )

    _object_terms: Set[str] = PrivateAttr(
        default={"http://semanticscience.org/resource/SIO_000370"}
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

    def add_to_regions(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        url: Optional[str] = None,
        accession_id: Optional[str] = None,
        name: Optional[str] = None,
        id: Optional[str] = None,
        **kwargs,
    ) -> Region:
        """
        This method adds an object of type 'Region' to attribute regions

        Args:
            id (str): Unique identifier of the 'Region' object. Defaults to 'None'.
            start (): Start position of the site.. Defaults to None
            end (): End position of the site.. Defaults to None
            url (): URI of the annotation.. Defaults to None
            accession_id (): Accession ID of the annotation.. Defaults to None
            name (): A name of a sequence feature, e.g. the name of a feature. Defaults to None
        """

        params = {
            "start": start,
            "end": end,
            "url": url,
            "accession_id": accession_id,
            "name": name,
        }

        if id is not None:
            params["id"] = id

        obj = Region(**params)

        self.regions.append(obj)

        return self.regions[-1]
