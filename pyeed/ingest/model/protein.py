from typing import Annotated, Any

from neo4j import AsyncDriver
from pydantic import Field, field_validator

from pyeed.ingest.model.annotation import Annotation

from .goannotation import GOAnnotation
from .pyeedbase import BaseNode, LabelProperty
from .reaction import Reaction
from .taxon import Taxon


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

    async def relate_to_reaction(
        self, driver: AsyncDriver, reaction: Reaction | list[Reaction]
    ) -> None:
        await self._relate(driver, "CATALYZES", reaction)

    async def relate_to_taxon(self, driver: AsyncDriver, taxon: Taxon | list[Taxon]) -> None:
        await self._relate(driver, "ORIGINATES_FROM", taxon)

    async def relate_to_go_annotation(
        self, driver: AsyncDriver, go_annotation: GOAnnotation | list[GOAnnotation]
    ) -> None:
        await self._relate(driver, "HAS_GO_ANNOTATION", go_annotation)

    async def relate_to_annotation(
        self, driver: AsyncDriver, annotation: Annotation | list[Annotation]
    ) -> None:
        await self._relate(driver, "HAS_ANNOTATION", annotation)
