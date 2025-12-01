from typing import Annotated, Any, Self

from pydantic import Field
from rdkit import Chem
from rdkit.Chem import inchi

from .pyeedbase import BaseNode, LabelProperty


class Molecule(BaseNode):
    """Chemical molecule information."""

    inchi_key: Annotated[str, LabelProperty(index=True)] = Field(
        ...,
        description="InChI key",
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

    @classmethod
    def from_smiles(
        cls, smiles: str, name: str | None = None, custom: dict[str, Any] | None = None
    ) -> Self:
        """Create a Molecule from a SMILES string.

        Args:
            smiles: SMILES string.
            name: Molecule name.
            custom: Custom data.

        Returns:
            Molecule object.

        Raises:
            ValueError: If the SMILES string cannot be parsed.
        """

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Could not parse SMILES: {smiles}")
        inchi_str = str(inchi.MolToInchi(mol))  # type: ignore
        inchikey = str(inchi.InchiToInchiKey(inchi_str))  # type: ignore
        return cls(
            inchi_key=inchikey,
            name=name,
            smiles=smiles,
            inchi=inchi_str,
            custom=custom or {},
        )
