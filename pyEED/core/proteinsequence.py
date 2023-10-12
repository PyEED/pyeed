import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from Bio import SeqIO, Entrez
from .equivalence import Equivalence
from .domain import Domain
from .annotation import Annotation
from .organism import Organism
from ..io_handler.sequence import _seqio_to_protein_sequence


@forge_signature
class ProteinSequence(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("proteinsequenceINDEX"),
        xml="@id",
    )

    name: str = Field(
        ...,
        description="Systematic name of the protein.",
    )

    amino_acid_sequence: str = Field(
        ...,
        description="The amino acid sequence of the protein sequence object.",
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

    organism: Organism = Field(
        ...,
        description="Corresponding organism",
    )

    domains: List[Domain] = Field(
        description="Domain specification",
        default_factory=ListPlus,
        multiple=True,
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

    annotations: List[Annotation] = Field(
        description="Position-wise annotation of the amino acid seqeunce",
        default_factory=ListPlus,
        multiple=True,
    )

    def add_to_domains(
        self,
        name: str,
        start_position: int,
        end_position: int,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Domain' to attribute domains

        Args:
            id (str): Unique identifier of the 'Domain' object. Defaults to 'None'.
            name (): Name of the annotated domain.
            start_position (): Position in the sequence where the domain starts.
            end_position (): Position in the sequence where the domain ends.
        """
        params = {
            "name": name,
            "start_position": start_position,
            "end_position": end_position,
        }
        if id is not None:
            params["id"] = id
        self.domains.append(Domain(**params))
        return self.domains[-1]

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

    def add_to_annotations(
        self,
        start_position: int,
        end_position: int,
        note: Optional[str] = None,
        name: Optional[str] = None,
        db_xref: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Annotation' to attribute annotations

        Args:
            id (str): Unique identifier of the 'Annotation' object. Defaults to 'None'.
            start_position (): Start position of the annotation. A single start position without an end corresponds to a single amino acid.
            end_position (): Optional end position if the annoation contains more than a single amino acid..
            note (): Function that is found in the annotated amino acid or. Defaults to None
            name (): Additional note for the annotation. Defaults to None
            db_xref (): Database cross reference. Defaults to None
        """
        params = {
            "start_position": start_position,
            "end_position": end_position,
            "note": note,
            "name": name,
            "db_xref": db_xref,
        }
        if id is not None:
            params["id"] = id
        self.annotations.append(Annotation(**params))
        return self.annotations[-1]

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
