from typing import Annotated

from pydantic import Field

from .pyeedbase import LabelProperty, PyeedBase


class Molecule(PyeedBase):
    """Chemical molecule information."""

    id: Annotated[str, LabelProperty(unique=True, index=True)] = Field(
        ...,
        description="Molecule identifier (ChEBI ID)",
    )
    name: str | None = Field(
        default=None,
        description="Molecule name",
    )
    smiles: str | None = Field(
        default=None,
        description="SMILES representation",
    )
    inchi: str | None = Field(
        default=None,
        description="InChI representation",
    )
