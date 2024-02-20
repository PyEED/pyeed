import sdRDM

from typing import Dict, List, Optional
from pydantic import PrivateAttr, model_validator
from uuid import uuid4
from pydantic_xml import attr, element, wrapped
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from .span import Span


@forge_signature
class AbstractRegion(
    sdRDM.DataModel,
    nsmap={
        "": "https://github.com/PyEED/pyeed@3b002efa6bf51e951767d8a7749ebad563897cb8#AbstractRegion"
    },
):
    """Annotation of a region within a sequence 🗺️"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    name: Optional[str] = element(
        description="Name of the annotation",
        default=None,
        tag="name",
        json_schema_extra=dict(),
    )

    spans: List[Span] = wrapped(
        "spans",
        element(
            description="Spans of the region. E.g. multiple exons of a gene",
            default_factory=ListPlus,
            tag="Span",
            json_schema_extra=dict(multiple=True),
        ),
    )

    note: Optional[str] = element(
        description="Information found in 'note' of an ncbi entry",
        default=None,
        tag="note",
        json_schema_extra=dict(),
    )

    cross_reference: Optional[str] = element(
        description="Database cross reference",
        default=None,
        tag="cross_reference",
        json_schema_extra=dict(),
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

    def add_to_spans(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        id: Optional[str] = None,
    ) -> Span:
        """
        This method adds an object of type 'Span' to attribute spans

        Args:
            id (str): Unique identifier of the 'Span' object. Defaults to 'None'.
            start (): Start position of the span of a region. Defaults to None
            end (): End position of the span of a region. Defaults to None
        """
        params = {"start": start, "end": end}
        if id is not None:
            params["id"] = id
        self.spans.append(Span(**params))
        return self.spans[-1]
