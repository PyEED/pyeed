import sdRDM

import time
import secrets
from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from requests import HTTPError
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from tqdm import tqdm
from .region import Region
from .site import Site
from .organism import Organism
from .nucleotidesequence import NucleotideSequence
from ..ncbi.seq_io import _seqio_to_protein_sequence, extract_nucleotide_seq


@forge_signature
class ProteinSequence(sdRDM.DataModel):
    """Description of a protein sequence. Additionally, the `ProteinSequence` contains annotations for sites and regions of the protein sequence alongside information on the organism. Furthermore, the `ProteinSequence` contains information on the coding sequence of the protein sequence, which allows later retrieval of the corresponding nucleotide sequence."""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("proteinsequenceINDEX"),
        xml="@id",
    )

    name: str = Field(
        ...,
        description="Name of the protein",
    )

    sequence: str = Field(
        ...,
        description="Amino acid sequence",
    )

    organism: Organism = Field(
        ...,
        description="Corresponding organism",
    )

    regions: List[Region] = Field(
        description="Domains of the protein",
        default_factory=ListPlus,
        multiple=True,
    )

    sites: List[Site] = Field(
        description="Annotations of different sites",
        default_factory=ListPlus,
        multiple=True,
    )

    coding_sequence: Optional[NucleotideSequence] = Field(
        default=NucleotideSequence(),
        description="Information about the coding sequence for the protein sequence",
    )

    ec_number: Optional[str] = Field(
        default=None,
        description="Enzyme Commission number",
    )

    mol_weight: Optional[float] = Field(
        default=None,
        description="Calculated molecular weight of the protein",
    )

    nr_id: Optional[str] = Field(
        default=None,
        description="Identifier for the NCBI NR database",
    )

    uniprot_id: Optional[str] = Field(
        default=None,
        description="Identifier for the UniProt database",
    )

    pdb_id: Optional[str] = Field(
        default=None,
        description="Identifier for the PDB database",
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

    def add_to_sites(
        self,
        name: Optional[str] = None,
        type: Optional[str] = None,
        positions: List[int] = ListPlus(),
        cross_reference: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Site' to attribute sites

        Args:
            id (str): Unique identifier of the 'Site' object. Defaults to 'None'.
            name (): Name of the site. Defaults to None
            type (): Type of the site. Defaults to None
            positions (): Positions of the site. Defaults to ListPlus()
            cross_reference (): Database cross reference. Defaults to None
        """
        params = {
            "name": name,
            "type": type,
            "positions": positions,
            "cross_reference": cross_reference,
        }
        if id is not None:
            params["id"] = id
        self.sites.append(Site(**params))
        return self.sites[-1]

    @classmethod
    def from_ncbi(cls, accession_id: str) -> "ProteinSequence":
        """
        This method creates a 'ProteinSequence' object from a given NCBI ID.

        Args:
            accession_id (str): NCBI accession ID of the protein sequence.

        Returns:
            ProteinSequence: 'ProteinSequence' object that corresponds to the given NCBI ID.
        """

        seq_record = cls._get_ncbi_entry(accession_id, "protein")
        return _seqio_to_protein_sequence(cls, seq_record)

    def pblast(self, n_hits: int) -> List["ProteinSequence"]:
        """Run protein blast for `ProteinSequence`.

        Args:
            n_hits (int): Number of hits to return.

        Returns:
            List[ProteinSequence]: List of 'ProteinSequence' objects that are the result of the blast search.
        """
        return self._pblast(self.sequence, n_hits)

    def _pblast(self, sequence: str, n_hits: int = None) -> List["ProteinSequence"]:
        print(f"Running pblast search for {self.name} from {self.organism.name}...")
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence, hitlist_size=n_hits)
        blast_record = NCBIXML.read(result_handle)

        accessions = self._get_accessions(blast_record)
        sequences = []
        bad_requests = []
        for acc in tqdm(accessions, desc="Fetching protein sequences"):
            try:
                sequences.append(ProteinSequence.from_ncbi(acc))
                time.sleep(1)

            except HTTPError:
                bad_requests.append(acc)

        if bad_requests:
            print(f"\nEncountered {len(bad_request)} bad HTTP requests, retrying...")
            for bad_request in tqdm(bad_requests, desc="Fetching protein sequences"):
                try:
                    sequences.append(ProteinSequence.from_ncbi(bad_request))
                    time.sleep(2)

                except HTTPError:
                    print(f"Failed to fetch {bad_request}")

        return sequences

    def _nblast(sequence: str, n_hits: int = None) -> List["ProteinSequence"]:
        result_handle = NCBIWWW.qblast("blastn", "nr", sequence, hitlist_size=n_hits)
        blast_record = NCBIXML.read(result_handle)

    @staticmethod
    def _get_accessions(blast_record: NCBIXML) -> List[str]:
        accessions = []
        for alignment in blast_record.alignments:
            accessions.append(alignment.accession)
        return accessions

    @property
    def nucleotide_seq(self):
        return self.coding_sequence.sequence
