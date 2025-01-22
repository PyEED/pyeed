from typing import List

from pydantic import BaseModel, Field


class Cluster(BaseModel):
    """
    Class for storing the cluster information
    """

    representative_id: str = Field(
        ...,
        description="ID of the representative sequence, representing sequences of a cluster.",
    )

    represented_ids: List[str] = Field(
        ...,
        description="List of IDs for sequences that are represented by the representative sequence in the cluster.",
    )


class Sequence(BaseModel):
    """
    Class for storing the sequence information
    """

    id: str = Field(..., description="ID of the sequence.")
    sequence: str = Field(..., description="Sequence of the sequence.")


class MultipleSequenceAlignment(BaseModel):
    """
    Class for storing the multiple sequence alignment with ids and sequences
    """

    sequences: List[Sequence] = Field(..., description="List of sequences.")

    def __str__(self) -> str:
        """Format alignment for display with sequences on separate lines."""
        if not self.sequences:
            return "Empty alignment"

        max_id_length = max(len(seq.id) for seq in self.sequences)
        return "\n".join(
            f"{seq.id:<{max_id_length}}  {seq.sequence}" for seq in self.sequences
        )

    def __repr__(self) -> str:
        return self.__str__()

    def __print__(self) -> None:
        print(self.__str__())
