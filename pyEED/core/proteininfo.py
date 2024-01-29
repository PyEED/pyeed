from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from Bio.Blast import NCBIWWW, NCBIXML
from pyEED.core.dnainfo import DNAInfo
from .proteinregion import ProteinRegion
from .abstractsequence import AbstractSequence
from .site import Site
from .citation import Citation
from .span import Span
from .proteinregiontype import ProteinRegionType
from .substrate import Substrate
from .dnaregion import DNARegion
from .proteinsitetype import ProteinSiteType
from ..ncbi.seq_io import _seqio_to_nucleotide_info, get_ncbi_entry, get_ncbi_entrys


@forge_signature
class ProteinInfo(AbstractSequence):
    """Description of a protein sequence. Additionally, the `ProteinSequence` contains annotations for sites and regions of the protein sequence alongside information on the organism. Furthermore, the `ProteinSequence` contains information on the coding sequence of the protein sequence, which allows later retrieval of the corresponding nucleotide sequence."""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("proteininfoINDEX"),
        xml="@id",
    )

    family_name: Optional[str] = Field(
        default=None,
        description="Family name of the protein",
    )

    regions: List[ProteinRegion] = Field(
        description="Domains of the protein",
        default_factory=ListPlus,
        multiple=True,
    )

    sites: List[Site] = Field(
        description="Annotations of different sites",
        default_factory=ListPlus,
        multiple=True,
    )

    coding_sequence_ref: Optional[DNARegion] = Field(
        description="Defines the coding sequence of the protein",
        default_factory=DNARegion,
    )

    ec_number: Optional[str] = Field(
        default=None,
        description="Enzyme Commission number",
    )

    mol_weight: Optional[float] = Field(
        default=None,
        description="Calculated molecular weight of the protein",
    )

    substrates: List[Substrate] = Field(
        description="Promiscuous substrates of the protein",
        default_factory=ListPlus,
        multiple=True,
    )

    def add_to_regions(
        self,
        type: Optional[ProteinRegionType] = None,
        name: Optional[str] = None,
        spans: List[Span] = ListPlus(),
        note: Optional[str] = None,
        cross_reference: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'ProteinRegion' to attribute regions

        Args:
            id (str): Unique identifier of the 'ProteinRegion' object. Defaults to 'None'.
            type (): Type of the region within the protein sequence. Defaults to None
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
        self.regions.append(ProteinRegion(**params))
        return self.regions[-1]

    def add_to_sites(
        self,
        name: Optional[str] = None,
        type: Optional[ProteinSiteType] = None,
        positions: List[int] = ListPlus(),
        cross_ref: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Site' to attribute sites

        Args:
            id (str): Unique identifier of the 'Site' object. Defaults to 'None'.
            name (): Name of the site. Defaults to None
            type (): Type of the site. Defaults to None
            positions (): Positions of the site. Defaults to ListPlus()
            cross_ref (): Database cross reference. Defaults to None
        """
        params = {
            "name": name,
            "type": type,
            "positions": positions,
            "cross_ref": cross_ref,
        }
        if id is not None:
            params["id"] = id
        self.sites.append(Site(**params))
        return self.sites[-1]

    def add_to_substrates(
        self,
        name: Optional[str] = None,
        inchi: Optional[str] = None,
        smiles: Optional[str] = None,
        chebi_id: Optional[str] = None,
        citation: Optional[Citation] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Substrate' to attribute substrates

        Args:
            id (str): Unique identifier of the 'Substrate' object. Defaults to 'None'.
            name (): Name of the substrate. Defaults to None
            inchi (): InChI code of the substrate. Defaults to None
            smiles (): SMILES code of the substrate. Defaults to None
            chebi_id (): ChEBI ID of the substrate. Defaults to None
            citation (): Citations of the substrate. Defaults to None
        """
        params = {
            "name": name,
            "inchi": inchi,
            "smiles": smiles,
            "chebi_id": chebi_id,
            "citation": citation,
        }
        if id is not None:
            params["id"] = id
        self.substrates.append(Substrate(**params))
        return self.substrates[-1]

    @classmethod
    def from_ncbi(cls, accession_id: str) -> "ProteinInfo":
        """
        This method creates a 'ProteinInfo' object from a given NCBI ID.

        Args:
            accession_id (str): NCBI accession ID of the protein sequence.

        Returns:
            ProteinInfo: 'ProteinInfo' object that corresponds to the given NCBI ID.
        """

        seq_record = get_ncbi_entry(accession_id, "protein")
        return _seqio_to_nucleotide_info(cls, seq_record)

    @classmethod
    def _from_seq_record(cls, seq_record) -> "ProteinInfo":
        return _seqio_to_nucleotide_info(cls, seq_record)

    def pblast(
        self,
        n_hits: int,
        e_value: float = 10.0,
        api_key: str = None,
        **kwargs,
    ) -> List["ProteinInfo"]:
        """Run protein blast for a `ProteinInfo`.
        Additional keyword arguments can be pass according to the blast [specifications](https://biopython.org/docs/1.75/api/Bio.Blast.NCBIWWW.html).

        Args:
            n_hits (int): Number of hits to return.
            e_value (float, optional): E-value threshold. Defaults to 10.0.
            api_key (str, optional): NCBI API key for sequence retrieval. Defaults to None.


        Returns:
            List[ProteinInfo]: List of 'ProteinInfo' objects that are the result of the blast search.
        """

        print(f"🏃🏼‍♀️ Running PBLAST")
        print(f"╭── protein name: {self.name}")
        print(f"├── accession: {self.source_id}")
        print(f"├── organism: {self.organism.name}")
        print(f"├── e-value: {e_value}")
        print(f"╰── max hits: {n_hits}")

        result_handle = NCBIWWW.qblast(
            "blastp",
            "nr",
            self.sequence,
            hitlist_size=n_hits,
            expect=e_value,
            **kwargs,
        )
        blast_record = NCBIXML.read(result_handle)

        accessions = self._get_accessions(blast_record)
        seq_records = get_ncbi_entrys(accessions, "protein", api_key=api_key)

        protein_infos = []
        for record in seq_records:
            protein_infos.append(self._from_seq_record(record))

        print("🎉 Done\n")
        return protein_infos

    def get_dna(self):
        if not self.coding_sequence_ref:
            return

        return DNAInfo.from_ncbi(self.coding_sequence_ref.id)

    def __str__(self):
        return f">{self.source_id}\n{self.sequence}"

    def _nblast(sequence: str, n_hits: int = None) -> List["ProteinInfo"]:
        result_handle = NCBIWWW.qblast("blastn", "nr", sequence, hitlist_size=n_hits)
        blast_record = NCBIXML.read(result_handle)

    @staticmethod
    def _get_accessions(blast_record: NCBIXML) -> List[str]:
        accessions = []
        for alignment in blast_record.alignments:
            accessions.append(alignment.accession)
        return accessions
