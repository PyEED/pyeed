import sdRDM

from typing import Optional, Union
from pydantic import PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from pydantic import Field
from typing import List
from typing import Optional

from .annotation import Annotation
from .domain import Domain
from .equivalence import Equivalence
from .organism import Organism


@forge_signature
class ProteinSequence(sdRDM.DataModel):
    id: str = Field(
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
        description="Identifier for the NCBI NR database",
        default=None,
    )

    uniprot_id: Optional[str] = Field(
        description="Identifier for the UniProt database",
        default=None,
    )

    pdb_id: List[str] = Field(
        description="Identifier for the PDB database",
        default_factory=ListPlus,
    )

    organism: Optional[Organism] = Field(
        description="Corresponding organism",
        default=None,
    )

    domain: List[Domain] = Field(
        description="Domain specification",
        default_factory=ListPlus,
    )

    reference_sequence: Optional[str] = Field(
        description="Identifier of the sequence used as reference",
        default=None,
    )

    equivalence: List[Equivalence] = Field(
        description="Positions where the given sequence is equivalent to the reference",
        default_factory=ListPlus,
    )

    annotation: List[Annotation] = Field(
        description="Position-wise annotation of the amino acid seqeunce",
        default_factory=ListPlus,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/PyEED/pyeed-data-model.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="047f17317fa860206980a47dc3790cbc3204f343"
    )

    def add_to_domain(
        self,
        name: str,
        start_position: int,
        end_position: int,
    ) -> None:
        """
        Adds an instance of 'Domain' to the attribute 'domain'.

        Args:
            name (str): Name of the annotated domain.
            start_position (int): Position in the sequence where the domain starts.
            end_position (int): Position in the sequence where the domain ends.
        """

        domain = [
            Domain(
                name=name,
                start_position=start_position,
                end_position=end_position,
            )
        ]

        self.domain = self.domain + domain

    def add_to_equivalence(
        self,
        reference_position: int,
        sequence_position: int,
    ) -> None:
        """
        Adds an instance of 'Equivalence' to the attribute 'equivalence'.

        Args:
            reference_position (int): Equivalent position in the reference sequence.
            sequence_position (int): Position that is equivalent to the reference sequence position that is also given.
        """

        equivalence = [
            Equivalence(
                reference_position=reference_position,
                sequence_position=sequence_position,
            )
        ]

        self.equivalence = self.equivalence + equivalence

    def add_to_annotation(
        self,
        start_position: int,
        function: str,
        end_position: Optional[int] = None,
    ) -> None:
        """
        Adds an instance of 'Annotation' to the attribute 'annotation'.

        Args:
            start_position (int): Start position of the annotation. A single start position without an end corresponds to a single amino acid.
            function (str): Function that is found in the annotated amino acid or sub-sequence.
            end_position (Optional[int]): Optional end position if the annoation contains more than a single amino acid. Defaults to None
        """

        annotation = [
            Annotation(
                start_position=start_position,
                function=function,
                end_position=end_position,
            )
        ]

        self.annotation = self.annotation + annotation
