from typing import Dict, List, Optional
from uuid import uuid4

import sdRDM
import validators
from lxml.etree import _Element
from pydantic import PrivateAttr, field_validator, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict


@forge_signature
class BlastData(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """"""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    identity: Optional[float] = element(
        description="Minimum identity to safe hits.",
        default=0.0,
        tag="identity",
        json_schema_extra=dict(),
    )

    evalue: Optional[float] = element(
        description="Expectation value (E) to safe hits.",
        default=10.0,
        tag="evalue",
        json_schema_extra=dict(),
    )

    n_hits: Optional[int] = element(
        description="Number of hits to return.",
        default=100,
        tag="n_hits",
        json_schema_extra=dict(),
    )

    substitution_matrix: Optional[str] = element(
        description="Substitution matrix to use.",
        default="'BLOSUM62'",
        tag="substitution_matrix",
        json_schema_extra=dict(),
    )

    word_size: Optional[int] = element(
        description="Word size of the initial match.",
        default=3,
        tag="word_size",
        json_schema_extra=dict(
            inclusivminimum=2,
            inclusivemaximum=7,
        ),
    )

    gap_open: Optional[float] = element(
        description="Gap open penalty.",
        default=11.0,
        tag="gap_open",
        json_schema_extra=dict(),
    )

    gap_extend: Optional[float] = element(
        description="Gap extend penalty.",
        default=1.0,
        tag="gap_extend",
        json_schema_extra=dict(),
    )

    threshold: Optional[float] = element(
        description="Minimum score to add a word to the BLAST lookup table.",
        default=11,
        tag="threshold",
        json_schema_extra=dict(),
    )

    db_name: Optional[str] = element(
        description="Name of the database to search.",
        default=None,
        tag="db_name",
        json_schema_extra=dict(),
    )

    annotations_: List[str] = element(
        tag="annotations_",
        alias="@type",
        description="Annotation of the given object.",
        default=[
            "BlastData",
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