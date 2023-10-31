import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .dnaregion import DNARegion
from .dnaregiontype import DNARegionType
from .organism import Organism


@forge_signature
class DNAInfo(sdRDM.DataModel):
    """Description of a nucleotide sequence 🧬"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("dnainfoINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the nucleotide sequence",
    )

    sequence: str = Field(
        ...,
        description="The nucleotide sequence coding for the protein sequence",
    )

    organism: Optional[Organism] = Field(
        default=None,
        description="Corresponding organism",
    )

    regions: List[DNARegion] = Field(
        description=(
            "Defines regions within the nucleotide sequence that code for the protein"
            " sequence"
        ),
        default_factory=ListPlus,
        multiple=True,
    )

    source_id: Optional[str] = Field(
        default=None,
        description="Identifier of the corresponding DNA sequence",
    )

    def add_to_regions(
        self,
        start: int,
        end: int,
        type: Optional[DNARegionType] = None,
        name: Optional[str] = None,
        note: Optional[str] = None,
        cross_reference: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'DNARegion' to attribute regions

        Args:
            id (str): Unique identifier of the 'DNARegion' object. Defaults to 'None'.
            start (): Start position of the annotation.
            end (): End position of the annotation.
            type (): Type of the region within the nucleotide sequence. Defaults to None
            name (): Name of the annotation. Defaults to None
            note (): Information found in 'note' of an ncbi protein sequence entry. Defaults to None
            cross_reference (): Database cross reference. Defaults to None
        """
        params = {
            "start": start,
            "end": end,
            "type": type,
            "name": name,
            "note": note,
            "cross_reference": cross_reference,
        }
        if id is not None:
            params["id"] = id
        self.regions.append(DNARegion(**params))
        return self.regions[-1]
