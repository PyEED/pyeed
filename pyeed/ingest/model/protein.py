from typing import Annotated, Any

from pydantic import Field, field_validator

from .pyeedbase import BaseNode, LabelProperty


class Protein(BaseNode):
    """Protein sequence and metadata."""

    id: Annotated[
        str,
        LabelProperty(index=True),
    ] = Field(
        ...,
        description="Protein identifier (UniProt accession)",
    )
    sequence: str = Field(
        ...,
        description="Amino acid sequence",
    )
    seq_length: int = Field(
        ...,
        description="Sequence length",
    )
    name: str | None = Field(
        None,
        description="Protein name",
    )
    mol_weight: float | None = Field(
        None,
        description="Molecular weight in Daltons",
    )
    ec_numbers: list[str] | None = Field(
        None,
        description="Enzyme Commission numbers associated with the protein",
    )

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("Sequence cannot be empty")
        return v.upper()

    @field_validator("seq_length")
    @classmethod
    def validate_seq_length(cls, v: int | None, info: Any) -> int | None:
        """Validate sequence length if provided, or auto-calculate from sequence."""
        sequence = info.data.get("sequence", "")
        if not sequence:
            return v

        actual_length = len(sequence)

        if v is None:
            # Auto-calculate length from sequence
            return actual_length
        elif v != actual_length:
            # Validate that provided length matches actual sequence length
            raise ValueError(
                f"Provided sequence length {v} does not match actual length {actual_length}"
            )
        elif v <= 0:
            raise ValueError("Sequence length must be positive")

        return v
