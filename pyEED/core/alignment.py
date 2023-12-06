import sdRDM

from typing import Optional, Union, List
from pydantic import PrivateAttr, Field, validator
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .span import Span
from .citation import Citation
from .author import Author
from .proteinregion import ProteinRegion
from .proteinregiontype import ProteinRegionType
from .organism import Organism
from .site import Site
from .standardnumbering import StandardNumbering
from .dnaregion import DNARegion
from .dnaregiontype import DNARegionType
from .substrate import Substrate
from .proteinsitetype import ProteinSiteType
from .proteininfo import ProteinInfo


@forge_signature
class Alignment(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("alignmentINDEX"),
        xml="@id",
    )

    reference_seq: Optional[ProteinInfo] = Field(
        default=None,
        description="Protein sequence used as reference",
        alias="reference",
    )

    query_seqs: List[ProteinInfo] = Field(
        description="Protein sequence used as query",
        default_factory=ListPlus,
        multiple=True,
    )

    method: Optional[str] = Field(
        default=None,
        description="Method used for the alignment",
    )

    consensus: Optional[str] = Field(
        default=None,
        description="Consensus sequence of the alignment",
    )

    score: Optional[float] = Field(
        default=None,
        description="Alignment score",
    )

    standard_numberings: List[StandardNumbering] = Field(
        description="Standard numbering of the aligned sequences",
        default_factory=ListPlus,
        multiple=True,
    )

    identity: Optional[float] = Field(
        default=None,
        description="Ration of identical residues in the alignment",
    )

    similarity: Optional[float] = Field(
        default=None,
        description="Ration of similar residues in the alignment",
    )

    gaps: Optional[int] = Field(
        default=None,
        description="Number of gaps in the alignment",
    )

    mismatches: Optional[int] = Field(
        default=None,
        description="Number of mismatches in the alignment",
    )

    def add_to_query_seqs(
        self,
        sequence: str,
        organism: Organism,
        source_id: Optional[str] = None,
        name: Optional[str] = None,
        regions: List[ProteinRegion] = ListPlus(),
        sites: List[Site] = ListPlus(),
        coding_sequence_ref: Optional[DNARegion] = None,
        ec_number: Optional[str] = None,
        mol_weight: Optional[float] = None,
        substrates: List[Substrate] = ListPlus(),
        citation: Optional[Citation] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'ProteinInfo' to attribute query_seqs

        Args:
            id (str): Unique identifier of the 'ProteinInfo' object. Defaults to 'None'.
            sequence (): Amino acid sequence.
            organism (): Corresponding organism.
            source_id (): Identifier of the protein sequence in the source database. Defaults to None
            name (): Name of the protein. Defaults to None
            regions (): Domains of the protein. Defaults to ListPlus()
            sites (): Annotations of different sites. Defaults to ListPlus()
            coding_sequence_ref (): Defines the coding sequence of the protein. Defaults to None
            ec_number (): Enzyme Commission number. Defaults to None
            mol_weight (): Calculated molecular weight of the protein. Defaults to None
            substrates (): Promiscuous substrates of the protein. Defaults to ListPlus()
            citation (): Publication on the protein. Defaults to None
        """
        params = {
            "sequence": sequence,
            "organism": organism,
            "source_id": source_id,
            "name": name,
            "regions": regions,
            "sites": sites,
            "coding_sequence_ref": coding_sequence_ref,
            "ec_number": ec_number,
            "mol_weight": mol_weight,
            "substrates": substrates,
            "citation": citation,
        }
        if id is not None:
            params["id"] = id
        self.query_seqs.append(ProteinInfo(**params))
        return self.query_seqs[-1]

    def add_to_standard_numberings(
        self,
        sequence_id: Optional[str] = None,
        numbering: List[str] = ListPlus(),
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'StandardNumbering' to attribute standard_numberings

        Args:
            id (str): Unique identifier of the 'StandardNumbering' object. Defaults to 'None'.
            sequence_id (): Identifier of the aligned sequence. Defaults to None
            numbering (): Standard numbering of the aligned sequence. Defaults to ListPlus()
        """
        params = {"sequence_id": sequence_id, "numbering": numbering}
        if id is not None:
            params["id"] = id
        self.standard_numberings.append(StandardNumbering(**params))
        return self.standard_numberings[-1]
