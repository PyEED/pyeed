import sdRDM

from typing import Dict, List, Optional
from pydantic import PrivateAttr, model_validator
from uuid import uuid4
from pydantic_xml import attr, element, wrapped
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict


@forge_signature
class Structure(
    sdRDM.DataModel,
    nsmap={
        "": "https://github.com/PyEED/pyeed@3b002efa6bf51e951767d8a7749ebad563897cb8#Structure"
    },
):
    """"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    pdb_id: Optional[str] = element(
        description="PDB ID of the structure",
        default=None,
        tag="pdb_id",
        json_schema_extra=dict(),
    )

    alphafold_id: Optional[str] = element(
        description="AlphaFold ID of the structure",
        default=None,
        tag="alphafold_id",
        json_schema_extra=dict(),
    )

    method: Optional[str] = element(
        description="Method used for structure determination",
        default=None,
        tag="method",
        json_schema_extra=dict(),
    )

    resolution: Optional[float] = element(
        description="Resolution of the structure in angstrom",
        default=None,
        tag="resolution",
        json_schema_extra=dict(),
    )

    chains: List[str] = wrapped(
        "chains",
        element(
            description="Chains of the structure",
            default_factory=ListPlus,
            tag="string",
            json_schema_extra=dict(multiple=True),
        ),
    )

    ligands: List[str] = wrapped(
        "ligands",
        element(
            description="Ligands of the structure",
            default_factory=ListPlus,
            tag="string",
            json_schema_extra=dict(multiple=True),
        ),
    )

    mutations: Optional[int] = element(
        description="Mutations of the structure",
        default=None,
        tag="mutations",
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
