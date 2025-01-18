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
