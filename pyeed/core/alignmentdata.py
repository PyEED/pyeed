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

from .sequence import Sequence


@forge_signature
class AlignmentData(
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

    consensus: Optional[str] = element(
        description="Consensus sequence of the alignment.",
        default=None,
        tag="consensus",
        json_schema_extra=dict(),
    )

    sequences: List[Sequence] = element(
        description="Sequences of the alignment.",
        default_factory=ListPlus,
        tag="sequences",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    aligned_sequences: List[Sequence] = element(
        description="Aligned sequences as a result of the alignment.",
        default_factory=ListPlus,
        tag="aligned_sequences",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    annotations_: List[str] = element(
        tag="annotations_",
        alias="@type",
        description="Annotation of the given object.",
        default=[
            "AlignmentData",
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

    def add_to_sequences(
        self,
        sequence_id: Optional[str] = None,
        sequence: Optional[str] = None,
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Sequence:
        """
        This method adds an object of type 'Sequence' to attribute sequences

        Args:
            id (str): Unique identifier of the 'Sequence' object. Defaults to 'None'.
            sequence_id (): Identifier of the sequence in the source database. Defaults to None
            sequence (): Molecular sequence.. Defaults to None
        """

        params = {
            "sequence_id": sequence_id,
            "sequence": sequence,
        }

        if id is not None:
            params["id"] = id

        obj = Sequence(**params)

        if annotation:
            obj.annotations_.append(annotation)

        self.sequences.append(obj)

        return self.sequences[-1]

    def add_to_aligned_sequences(
        self,
        sequence_id: Optional[str] = None,
        sequence: Optional[str] = None,
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Sequence:
        """
        This method adds an object of type 'Sequence' to attribute aligned_sequences

        Args:
            id (str): Unique identifier of the 'Sequence' object. Defaults to 'None'.
            sequence_id (): Identifier of the sequence in the source database. Defaults to None
            sequence (): Molecular sequence.. Defaults to None
        """

        params = {
            "sequence_id": sequence_id,
            "sequence": sequence,
        }

        if id is not None:
            params["id"] = id

        obj = Sequence(**params)

        if annotation:
            obj.annotations_.append(annotation)

        self.aligned_sequences.append(obj)

        return self.aligned_sequences[-1]
