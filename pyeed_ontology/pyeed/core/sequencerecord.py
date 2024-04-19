import sdRDM

import validators
from typing import Dict, List, Optional
from pydantic import PrivateAttr, field_validator, model_validator
from uuid import uuid4
from pydantic_xml import attr, element
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from .organism import Organism


@forge_signature
class SequenceRecord(sdRDM.DataModel, search_mode="unordered"):
    """A molecular sequence and associated annotation data."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    uri: Optional[str] = element(
        description="URI of the sequence.",
        default=None,
        tag="uri",
        json_schema_extra=dict(term="http://edamontology.org/data_0848"),
    )

    accession_id: Optional[str] = element(
        description="Accession ID of the sequence.",
        default=None,
        tag="accession_id",
        json_schema_extra=dict(term="http://edamontology.org/data_2091"),
    )

    name: Optional[str] = element(
        description="Arbtrary name of the sequence.",
        default=None,
        tag="name",
        json_schema_extra=dict(term="http://edamontology.org/data_1009"),
    )

    organism: Optional[Organism] = element(
        description="The organism from which the sequence was obtained.",
        default=None,
        tag="organism",
        json_schema_extra=dict(),
    )

    annotations_: List[str] = element(
        tag="annotations_",
        alias="@type",
        description="Annotation of the given object.",
        default=["SequenceRecord", "http://edamontology.org/data_0848"],
    )
    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="b6c966bc9237ce7813e3e8dd6e17f57734d84693"
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

    @field_validator("annotations_")
    @classmethod
    def _validate_annotation(cls, annotations):
        """Check if the annotation that has been set is a valid URL."""
        for annotation in annotations:
            if not validators.url(annotation) and annotation != cls.__name__:
                raise ValueError(f"Invalid Annotation URL: {annotation}")
        return annotations

    @classmethod
    def _validate_annotation(cls, annotations):
        """Check if the annotation that has been set is a valid URL."""

        for annotation in annotations:
            if not validators.url(annotation) and annotation != cls.__name__:
                raise ValueError(f"Invalid Annotation URL: {annotation}")

        return annotations
