from abc import ABC, abstractmethod
from typing import List

from pydantic import BaseModel, Field


class AbstractAligner(BaseModel, ABC):

    sequences: List[str] = Field(
        ..., description="List of sequence strings to be aligned"
    )

    def __init__(self, **kwargs):
        BaseModel.__init__(self, **kwargs)
        super().__init__(**kwargs)

    @abstractmethod
    def align(self):
        """
        Abstract method for aligning sequences.
        """
        pass
