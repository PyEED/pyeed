from typing import Annotated

from pydantic import Field

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


from rdkit import Chem
from rdkit.Chem import inchi


def smiles_to_inchikey(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    inchi_str = inchi.MolToInchi(mol)  # uses IUPAC InChI under the hood
    inchikey = inchi.InchiToInchiKey(inchi_str)
    return inchikey


smiles = "C1CCCCC1"
inchikey = smiles_to_inchikey(smiles)
print(inchikey)
