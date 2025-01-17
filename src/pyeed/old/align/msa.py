from Bio.Align import MultipleSeqAlignment
from pydantic import BaseModel, Field
from pyeed.core.alignmentresult import AlignmentResult
from pyeed.core.proteinrecord import ProteinRecord
from pyeed.core.sequence import Sequence
from pyeed.tools.clustalo import ClustalOmega
from rich.console import Console
from rich.status import Status


class MSA(BaseModel):
    sequences: list[Sequence | ProteinRecord] = Field(
        ...,
        description="List of sequences to be aligned.",
        min_length=2,
    )

    def clustalo(self) -> AlignmentResult:
        """Aligns the sequences using Clustal Omega.

        Returns:
            AlignmentResult: The alignment result.
        """

        assert (
            len(self.sequences) >= 2
        ), "At least two sequences are required for alignment."

        multifasta = [seq.fasta_string() for seq in self.sequences]

        with Status("Running ClustalOmega...", console=Console(force_terminal=False)):
            alignment = ClustalOmega().align(multifasta)
        print("✅ Alignment completed")

        return self._map_alignment(alignment)

    def _map_alignment(self, alignment: MultipleSeqAlignment) -> AlignmentResult:
        """Maps the alignment Biopython result to the AlignmentResult object.

        Args:
            alignment (MultipleSeqAlignment): Biopython result object.

        Returns:
            AlignmentResult: Result of the alignment.
        """

        result = AlignmentResult()
        [
            result.add_to_sequences(
                id=seq.id, sequence=seq.sequence, sequence_id=seq.id
            )
            for seq in self.sequences
        ]
        [
            result.add_to_aligned_sequences(
                id=seq.id, sequence=str(seq.seq), sequence_id=seq.id
            )
            for seq in alignment
        ]

        return result


if __name__ == "__main__":
    from pyeed.core.proteinrecord import ProteinRecord
    from pyeed.core.sequence import Sequence

    protein1 = ProteinRecord(
        id="seq1",
        sequence="MTHKLLLTLLFTLLFSSAYSRG",
    )
    protein2 = ProteinRecord(
        id="seq2",
        sequence="MTHKILLLTLLFTLLFSDSSAYSRG",
    )
    protein3 = ProteinRecord(
        id="seq3",
        sequence="MTHKILLLTLLFTLLFSSCYSRG",
    )
    protein4 = Sequence(
        id="seq4",
        sequence="AAAAAAA",
    )

    sequences = [protein1, protein2, protein3, protein4]

    msa = MSA(sequences=sequences)
    alignment = msa.clustalo()
    print(alignment)
