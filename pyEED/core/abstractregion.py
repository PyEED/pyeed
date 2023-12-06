import sdRDM

from typing import Optional, Union, List
from pydantic import PrivateAttr, Field, validator
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .span import Span


@forge_signature
class AbstractRegion(sdRDM.DataModel):
    """Annotation of a region within a sequence 🗺️"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("abstractregionINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the annotation",
    )

    spans: List[Span] = Field(
        description="Spans of the region. E.g. multiple exons of a gene",
        default_factory=ListPlus,
        multiple=True,
    )

    note: Optional[str] = Field(
        default=None,
        description="Information found in 'note' of an ncbi entry",
    )

    cross_reference: Optional[str] = Field(
        default=None,
        description="Database cross reference",
    )

    def add_to_spans(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        id: Optional[str] = None,
    ) -> None:
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
