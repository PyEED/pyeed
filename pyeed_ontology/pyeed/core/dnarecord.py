import validators

from typing import Dict, List, Optional
from uuid import uuid4
from pydantic import PrivateAttr, field_validator, model_validator
from pydantic_xml import attr, element
from lxml.etree import _Element

from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict


from .region import Region
from .site import Site
from .sequencerecord import SequenceRecord


@forge_signature
class DNARecord(
    SequenceRecord,
    search_mode="unordered",
):
    """Description of a nucleotide sequence 🧬."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    sequence: str = element(
        description="Nucleotide sequence.",
        tag="sequence",
        json_schema_extra=dict(
            term="http://edamontology.org/data_3494",
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

    sites: List[Site] = element(
        description=(
            "Defines sites within the nucleotide sequence that code for the protein"
            " sequence."
        ),
        default_factory=ListPlus,
        tag="sites",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    gc_content: Optional[float] = element(
        description="GC content of the sequence.",
        default=None,
        tag="gc_content",
        json_schema_extra=dict(),
    )

    annotations_: List[str] = element(
        tag="annotations_",
        alias="@type",
        description="Annotation of the given object.",
        default=[
            "DNARecord",
        ],
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="f113b5b736593e06d2e2ded44e4a2c83052e2fbc"
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

    @field_validator("annotations_")
    @classmethod
    def _validate_annotation(cls, annotations):
        """Check if the annotation that has been set is a valid URL."""

        for annotation in annotations:
            if not validators.url(annotation) and annotation != cls.__name__:
                raise ValueError(f"Invalid Annotation URL: {annotation}")

        return annotations

    def add_to_regions(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Region:
        """
        This method adds an object of type 'Region' to attribute regions

        Args:
            id (str): Unique identifier of the 'Region' object. Defaults to 'None'.
            start (): Start position of the site.. Defaults to None
            end (): End position of the site.. Defaults to None
        """

        params = {
            "start": start,
            "end": end,
        }

        if id is not None:
            params["id"] = id

        obj = Region(**params)

        if annotation:
            obj.annotations_.append(annotation)

        self.regions.append(obj)

        return self.regions[-1]

    def add_to_sites(
        self,
        positions: List[int] = ListPlus(),
        uri: Optional[str] = None,
        accession_id: Optional[str] = None,
        name: Optional[str] = None,
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Site:
        """
        This method adds an object of type 'Site' to attribute sites

        Args:
            id (str): Unique identifier of the 'Site' object. Defaults to 'None'.
            positions (): Position of the site(s) within the sequence.. Defaults to ListPlus()
            uri (): URI of the annotation.. Defaults to None
            accession_id (): Accession ID of the annotation.. Defaults to None
            name (): A name of a sequence feature, e.g. the name of a feature. Defaults to None
        """

        params = {
            "positions": positions,
            "uri": uri,
            "accession_id": accession_id,
            "name": name,
        }

        if id is not None:
            params["id"] = id

        obj = Site(**params)

        if annotation:
            obj.annotations_.append(annotation)

        self.sites.append(obj)

        return self.sites[-1]
