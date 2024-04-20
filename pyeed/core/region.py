from typing import Dict, List, Optional
from uuid import uuid4

import sdRDM
import validators
from lxml.etree import _Element
from pydantic import PrivateAttr, field_validator, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict


class Region(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """Regional annotation of a feature within a sequence. The direction of the region is defined by the start and end positions."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    start: Optional[int] = element(
        description="Start position of the site.",
        default=None,
        tag="start",
        json_schema_extra=dict(),
    )

    end: Optional[int] = element(
        description="End position of the site.",
        default=None,
        tag="end",
        json_schema_extra=dict(),
    )

    annotations_: List[str] = element(
        tag="annotations_",
        alias="@type",
        description="Annotation of the given object.",
        default=[
            "Region",
        ],
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
