import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .region import Region


@forge_signature
class NucleotideSequence(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("nucleotidesequenceINDEX"),
        xml="@id",
    )

    regions: List[Region] = Field(
        description=(
            "Defines regions within the nucleotide sequence that code for the protein"
            " sequence"
        ),
        default_factory=ListPlus,
        multiple=True,
    )

    molecule_type: Optional[str] = Field(
        default=None,
        description="Type of the sequence",
    )

    protein_id: Optional[str] = Field(
        default=None,
        description="Identifier of the corresponding protein sequence",
    )

    gene_id: Optional[str] = Field(
        default=None,
        description="Identifier of the corresponding gene",
    )

    sequence: Optional[str] = Field(
        default=None,
        description="The nucleotide sequence coding for the protein sequence",
    )

    def add_to_regions(
        self,
        start: int,
        end: int,
        note: Optional[str] = None,
        name: Optional[str] = None,
        cross_reference: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Region' to attribute regions

        Args:
            id (str): Unique identifier of the 'Region' object. Defaults to 'None'.
            start (): Start position of the annotation. A single start position without an end corresponds to a single amino acid.
            end (): Optional end position if the annotation contains more than a single amino acid.
            note (): Information found in 'note' of an ncbi protein sequence entry. Defaults to None
            name (): Name of the annotation. Defaults to None
            cross_reference (): Database cross reference. Defaults to None
        """
        params = {
            "start": start,
            "end": end,
            "note": note,
            "name": name,
            "cross_reference": cross_reference,
        }
        if id is not None:
            params["id"] = id
        self.regions.append(Region(**params))
        return self.regions[-1]
