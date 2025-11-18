from typing import Annotated, Any

from pydantic import Field, field_validator

from .pyeedbase import LabelProperty, PyeedBase


class Protein(PyeedBase):
    """Protein sequence and metadata."""

    id: Annotated[
        str,
        LabelProperty(unique=True, index=True),
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
    taxon_ids: list[str] = Field(
        default_factory=list,
        description="Taxonomy IDs of the organism the protein originates from",
    )
    structure_ids: list[str] = Field(
        default_factory=list,
        description="Structure identifiers",
    )
    go_ids: list[str] = Field(
        default_factory=list,
        description="Gene Ontology identifiers",
    )
    reaction_ids: list[str] = Field(
        default_factory=list,
        description="RHEA reaction identifiers",
    )
    annotation_ids: list[str] = Field(
        default_factory=list,
        description="Sequence annotation identifiers",
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
