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
from .doi import DOI


@forge_signature
class Annotation(sdRDM.DataModel, search_mode="unordered"):
    """"""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    accession_id: Optional[str] = element(
        description="Accession ID of the annotation.",
        default=None,
        tag="accession_id",
        json_schema_extra=dict(term="http://edamontology.org/data_3034"),
    )

    name: Optional[str] = element(
        description="A name of a sequence feature, e.g. the name of a feature",
        default=None,
        tag="name",
        json_schema_extra=dict(),
    )

    publications: List[DOI] = element(
        description="DOI of the publication(s) referenced in the annotation.",
        default_factory=ListPlus,
        tag="publications",
        json_schema_extra=dict(multiple=True),
    )

    annotations_: List[str] = element(
        tag="annotations_",
        alias="@type",
        description="Annotation of the given object.",
        default=["Annotation"],
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

    def add_to_publications(
        self,
        doi: Optional[str] = None,
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> DOI:
        """
        This method adds an object of type 'DOI' to attribute publications

        Args:
            id (str): Unique identifier of the 'DOI' object. Defaults to 'None'.
            doi (): Digital Object Identifier (DOI) of a publication.. Defaults to None
        """
        params = {"doi": doi}
        if id is not None:
            params["id"] = id
        obj = DOI(**params)
        if annotation:
            obj.annotations_.append(annotation)
        self.publications.append(obj)
        return self.publications[-1]

    @classmethod
    def _validate_annotation(cls, annotations):
        """Check if the annotation that has been set is a valid URL."""

        for annotation in annotations:
            if not validators.url(annotation) and annotation != cls.__name__:
                raise ValueError(f"Invalid Annotation URL: {annotation}")

        return annotations

    def add_to_publications(
        self,
        doi: Optional[str] = None,
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> DOI:
        """
        This method adds an object of type 'DOI' to attribute publications

        Args:
            id (str): Unique identifier of the 'DOI' object. Defaults to 'None'.
            doi (): Digital Object Identifier (DOI) of a publication.. Defaults to None
        """

        params = {
            "doi": doi,
        }

        if id is not None:
            params["id"] = id

        obj = DOI(**params)

        if annotation:
            obj.annotations_.append(annotation)

        self.publications.append(obj)

        return self.publications[-1]
