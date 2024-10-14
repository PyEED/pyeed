from typing import Dict, List, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from .organism import Organism
from .region import Region
from .regionset import RegionSet
from .site import Site


class SequenceRecord(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """A molecular sequence and associated annotation data."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    id: Optional[str] = element(
        description="Unique identifier of the sequence.",
        default=None,
        tag="id",
        json_schema_extra=dict(
            term="http://semanticscience.org/resource/SIO_000729",
        ),
    )

    name: Optional[str] = element(
        description="Arbitrary name of the sequence.",
        default=None,
        tag="name",
        json_schema_extra=dict(
            term="http://semanticscience.org/resource/SIO_000116",
        ),
    )

    organism: Optional[Organism] = element(
        description="The organism from which the sequence was obtained.",
        default=None,
        tag="organism",
        json_schema_extra=dict(
            term="http://semanticscience.org/resource/SIO_010000",
        ),
    )

    sequence: str = element(
        description="The letter sequence of the macromolecule.",
        tag="sequence",
        json_schema_extra=dict(
            term="http://semanticscience.org/resource/SIO_000030",
        ),
    )

    seq_length: Optional[int] = element(
        description="Length of the sequence.",
        default=None,
        tag="seq_length",
        json_schema_extra=dict(
            term="http://semanticscience.org/resource/SIO_000041",
        ),
    )

    sites: List[Site] = element(
        description="Defines sites within the nucleotide sequence.",
        default_factory=ListPlus,
        tag="sites",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    regions: List[Region] = element(
        description="Defines regions within the nucleotide sequence.",
        default_factory=ListPlus,
        tag="regions",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    region_sets: List[RegionSet] = element(
        description=(
            "Multiple regions forming a higher order structure or feature of a"
            " sequence."
        ),
        default_factory=ListPlus,
        tag="region_sets",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="b926bfec3aa1ec45a5614cf6ac4a546252dd384c"
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

    def add_to_sites(
        self,
        positions: List[int] = ListPlus(),
        url: Optional[str] = None,
        accession_id: Optional[str] = None,
        name: Optional[str] = None,
        id: Optional[str] = None,
        **kwargs,
    ) -> Site:
        """
        This method adds an object of type 'Site' to attribute sites

        Args:
            id (str): Unique identifier of the 'Site' object. Defaults to 'None'.
            positions (): Position of the site(s) within the sequence.. Defaults to ListPlus()
            url (): URI of the annotation.. Defaults to None
            accession_id (): Accession ID of the annotation.. Defaults to None
            name (): A name of a sequence feature, e.g. the name of a feature. Defaults to None
        """

        params = {
            "positions": positions,
            "url": url,
            "accession_id": accession_id,
            "name": name,
        }

        if id is not None:
            params["id"] = id

        obj = Site(**params)

        self.sites.append(obj)

        return self.sites[-1]

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

    def add_to_region_sets(
        self,
        regions: List[Region] = ListPlus(),
        id: Optional[str] = None,
        **kwargs,
    ) -> RegionSet:
        """
        This method adds an object of type 'RegionSet' to attribute region_sets

        Args:
            id (str): Unique identifier of the 'RegionSet' object. Defaults to 'None'.
            regions (): Regions of the cluster.. Defaults to ListPlus()
        """

        params = {
            "regions": regions,
        }

        if id is not None:
            params["id"] = id

        obj = RegionSet(**params)

        self.region_sets.append(obj)

        return self.region_sets[-1]

    def fasta_string(self):
        """Return the sequence as a fasta formatted string."""

        return f">{self.id}\n{self.sequence}\n"

    def to_fasta(self, path):
        """Write the sequence to a FASTA file."""

        with open(path, "w") as f:
            f.write(self.fasta_string())
