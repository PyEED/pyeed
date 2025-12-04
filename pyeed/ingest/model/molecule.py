from typing import Annotated

from pydantic import Field

from .pyeedbase import BaseNode, LabelProperty


class Molecule(BaseNode):
    """Chemical molecule information."""

    id: Annotated[str, LabelProperty(index=True)] = Field(
        ...,
        description="ChEBI ID",
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

    # @classmethod
    # def from_smiles(
    #     cls,
    #     smiles: str,
    #     name: str | None = None,
    # ) -> Self:
    #     """Create a Molecule from a SMILES string.

    #     Args:
    #         smiles: SMILES string.
    #         name: Molecule name.
    #         custom: Custom data.

    #     Returns:
    #         Molecule object.

    #     Raises:
    #         ValueError: If the SMILES string cannot be parsed.
    #     """

    #     mol = Chem.MolFromSmiles(smiles)
    #     if mol is None:
    #         raise ValueError(f"Could not parse SMILES: {smiles}")
    #     inchi_str = str(inchi.MolToInchi(mol))  # type: ignore
    #     inchikey = str(inchi.InchiToInchiKey(inchi_str))  # type: ignore
    #     return cls(
    #         inchi_key=inchikey,
    #         name=name,
    #         smiles=smiles,
    #         inchi=inchi_str,
    #     )
