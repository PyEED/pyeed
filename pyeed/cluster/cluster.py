from pydantic import BaseModel, Field
from pyeed.core import AbstractSequence


class Cluster(BaseModel):
    """Data model for a cluster of sequences."""

    name: str = Field(
        description="The name of the cluster",
        default=None,
    )
    representative: AbstractSequence = Field(
        description="The representative sequence of the cluster",
        default=None,
    )
    members: list[AbstractSequence] = Field(
        description="The sequences in the cluster",
        default=None,
    )
