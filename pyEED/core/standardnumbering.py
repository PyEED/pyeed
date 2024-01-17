import sdRDM

from typing import Optional, Union, List
from pydantic import Field, field_validator
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .abstractsequence import AbstractSequence


@forge_signature
class StandardNumbering(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("standardnumberingINDEX"),
        xml="@id",
    )

    sequence_id: Union[AbstractSequence, str, None] = Field(
        default=None,
        reference="AbstractSequence.source_id",
        description="Identifier of the aligned sequence",
    )

    numbering: List[str] = Field(
        description="Standard numbering of the aligned sequence",
        default_factory=ListPlus,
        multiple=True,
    )

    @field_validator("sequence_id")
    def get_sequence_id_reference(cls, value):
        """Extracts the ID from a given object to create a reference"""
        from .abstractsequence import AbstractSequence

        if isinstance(value, AbstractSequence):
            return value.source_id
        elif isinstance(value, str):
            return value
        elif value is None:
            return value
        else:
            raise TypeError(
                f"Expected types [AbstractSequence, str] got '{type(value).__name__}'"
                " instead."
            )
