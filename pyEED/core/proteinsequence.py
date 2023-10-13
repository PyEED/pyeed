import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from Bio import SeqIO, Entrez
from .organism import Organism
from .equivalence import Equivalence
from .dnasequence import DNASequence
from .annotation import Annotation
from ..io_handler.sequence import _seqio_to_protein_sequence


@forge_signature
class ProteinSequence(sdRDM.DataModel):
    """Description of a protein sequence and its annotations."""

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

    regions: List[Annotation] = Field(
        description="Domains of the protein",
        default_factory=ListPlus,
        multiple=True,
    )

    sites: List[Annotation] = Field(
        description="Annotations of different sites",
        default_factory=ListPlus,
        multiple=True,
    )

    cds: Optional[DNASequence] = Field(
        default=None,
        description="Corresponding DNA coding sequence",
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

    reference_sequence: Optional[str] = Field(
        default=None,
        description="Identifier of the sequence used as reference",
    )

    equivalence: List[Equivalence] = Field(
        description="Positions where the given sequence is equivalent to the reference",
        default_factory=ListPlus,
        multiple=True,
    )
    __repo__: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed.git")
    __commit__: Optional[str] = PrivateAttr(
        default="5592833c9e695a2fe5a0070d25ddb41fdeab9f05"
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
        This method adds an object of type 'Annotation' to attribute regions

        Args:
            id (str): Unique identifier of the 'Annotation' object. Defaults to 'None'.
            start (): Start position of the annotation. A single start position without an end corresponds to a single amino acid.
            end (): Optional end position if the annoation contains more than a single amino acid.
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
        self.regions.append(Annotation(**params))
        return self.regions[-1]

    def add_to_sites(
        self,
        start: int,
        end: int,
        note: Optional[str] = None,
        name: Optional[str] = None,
        cross_reference: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Annotation' to attribute sites

        Args:
            id (str): Unique identifier of the 'Annotation' object. Defaults to 'None'.
            start (): Start position of the annotation. A single start position without an end corresponds to a single amino acid.
            end (): Optional end position if the annoation contains more than a single amino acid.
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
        self.sites.append(Annotation(**params))
        return self.sites[-1]

    def add_to_equivalence(
        self, reference_position: int, sequence_position: int, id: Optional[str] = None
    ) -> None:
        """
        This method adds an object of type 'Equivalence' to attribute equivalence

        Args:
            id (str): Unique identifier of the 'Equivalence' object. Defaults to 'None'.
            reference_position (): Equivalent position in the reference sequence.
            sequence_position (): Position that is equivalent to the reference sequence position that is also given.
        """
        params = {
            "reference_position": reference_position,
            "sequence_position": sequence_position,
        }
        if id is not None:
            params["id"] = id
        self.equivalence.append(Equivalence(**params))
        return self.equivalence[-1]

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

    def _get_ncbi_entry(
        accession_id: str, database: str, email: str = "exon@got-spliced.com"
    ) -> SeqIO.SeqRecord:
        Entrez.email = email

        databases = {"nucleotide", "protein"}
        if database not in databases:
            raise ValueError(f"database must be one of {databases}")

        handle = Entrez.efetch(
            db=database, id=accession_id, rettype="gb", retmode="text"
        )
        seq_record = SeqIO.read(handle, "genbank")

        handle.close()
        return seq_record
