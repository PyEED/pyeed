from typing import Dict, List, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from .sequence import Sequence


class Cluster(
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

    name: Optional[str] = element(
        description="Name of the cluster.",
        default=None,
        tag="name",
        json_schema_extra=dict(),
    )

    representative: Optional[Sequence] = element(
        description="Identifier of the representative sequence of the cluster.",
        default_factory=Sequence,
        tag="representative",
        json_schema_extra=dict(),
    )

    members: List[Sequence] = element(
        description="Sequences of the cluster.",
        default_factory=ListPlus,
        tag="members",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="58f6e926b555159b778f5248737b8d20ea09fca0"
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

    def add_to_members(
        self,
        sequence_id: Optional[str] = None,
        sequence: Optional[str] = None,
        id: Optional[str] = None,
        **kwargs,
    ) -> Sequence:
        """
        This method adds an object of type 'Sequence' to attribute members

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

        self.members.append(obj)

        return self.members[-1]
