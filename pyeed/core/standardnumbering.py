from typing import Dict, List, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from .numberedsequence import NumberedSequence


class StandardNumbering(
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

    reference_id: Optional[str] = element(
        description="Standard numbering of the reference sequence",
        default=None,
        tag="reference_id",
        json_schema_extra=dict(),
    )

    numberd_sequences: List[NumberedSequence] = element(
        description="Numbered sequence of the aligned sequence",
        default_factory=ListPlus,
        tag="numberd_sequences",
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

    def add_to_numberd_sequences(
        self,
        numbered_id: Optional[str] = None,
        numbering: List[str] = ListPlus(),
        id: Optional[str] = None,
        **kwargs,
    ) -> NumberedSequence:
        """
        This method adds an object of type 'NumberedSequence' to attribute numberd_sequences

        Args:
            id (str): Unique identifier of the 'NumberedSequence' object. Defaults to 'None'.
            numbered_id (): Identifier of the numbered sequence. Defaults to None
            numbering (): Standard numbering of the aligned sequence. Defaults to ListPlus()
        """

        params = {
            "numbered_id": numbered_id,
            "numbering": numbering,
        }

        if id is not None:
            params["id"] = id

        obj = NumberedSequence(**params)

        self.numberd_sequences.append(obj)

        return self.numberd_sequences[-1]
