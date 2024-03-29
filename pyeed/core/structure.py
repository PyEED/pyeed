import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Structure(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("structureINDEX"),
        xml="@id",
    )

    pdb_id: Optional[str] = Field(
        default=None,
        description="PDB ID of the structure",
    )

    alphafold_id: Optional[str] = Field(
        default=None,
        description="AlphaFold ID of the structure",
    )

    method: Optional[str] = Field(
        default=None,
        description="Method used for structure determination",
    )

    resolution: Optional[float] = Field(
        default=None,
        description="Resolution of the structure in angstrom",
    )

    chains: List[str] = Field(
        description="Chains of the structure",
        default_factory=ListPlus,
        multiple=True,
    )

    ligands: List[str] = Field(
        description="Ligands of the structure",
        default_factory=ListPlus,
        multiple=True,
    )

    mutations: Optional[int] = Field(
        default=None,
        description="Mutations of the structure",
    )
