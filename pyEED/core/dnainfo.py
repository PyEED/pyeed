import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .dnaregion import DNARegion
from .organism import Organism
from .span import Span
from .dnaregiontype import DNARegionType
from ..ncbi.seq_io import get_ncbi_entry, _seqio_to_dna_info


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
        type: Optional[DNARegionType] = None,
        name: Optional[str] = None,
        spans: List[Span] = ListPlus(),
        note: Optional[str] = None,
        cross_reference: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'DNARegion' to attribute regions

        Args:
            id (str): Unique identifier of the 'DNARegion' object. Defaults to 'None'.
            type (): Type of the region within the nucleotide sequence. Defaults to None
            name (): Name of the annotation. Defaults to None
            spans (): Spans of the region. E.g. multiple exons of a gene. Defaults to ListPlus()
            note (): Information found in 'note' of an ncbi entry. Defaults to None
            cross_reference (): Database cross reference. Defaults to None
        """
        params = {
            "type": type,
            "name": name,
            "spans": spans,
            "note": note,
            "cross_reference": cross_reference,
        }
        if id is not None:
            params["id"] = id
        self.regions.append(DNARegion(**params))
        return self.regions[-1]

    @classmethod
    def from_ncbi(cls, accession_id: str) -> "DNAInfo":
        seq_record = get_ncbi_entry(accession_id=accession_id, database="nucleotide")
        return _seqio_to_dna_info(cls, seq_record)

    # def extract_nucleotide_seq(self, dna_region: DNARegion):
    #     """Handel nucleotide SeqIO entry and map it to `NucleotideSequence`"""

    #     for feature in entry.features:
    #         if feature.type == "CDS":
    #             if isinstance(feature.location, CompoundLocation):
    #                 locations = set()
    #                 parts = feature.location.parts
    #                 for part in parts:
    #                     # TODO: investigate reason for +1
    #                     locations.add(int(part.start) + 1)
    #                     locations.add(int(part.end))

    #             if isinstance(feature.location, FeatureLocation):
    #                 locations = {int(feature.location.start) + 1, int(feature.location.end)}

    #             if feature_regions == locations:
    #                 nucleotide_sequence.sequence = str(feature.location.extract(entry.seq))
