from __future__ import annotations

import strawberry


@strawberry.input(description="Filter proteins by any field")
class ProteinFilter:
    """Filter input for querying proteins by any field.

    All fields are optional. Multiple fields are combined with AND.
    For range queries (e.g., seqLength), use min/max variants.
    """

    id: str | None = None
    name: str | None = None
    sequence: str | None = None
    seq_length: int | None = None
    seq_length_min: int | None = strawberry.field(
        default=None, description="Minimum sequence length (inclusive)"
    )
    seq_length_max: int | None = strawberry.field(
        default=None, description="Maximum sequence length (inclusive)"
    )
    mol_weight: float | None = None
    mol_weight_min: float | None = strawberry.field(
        default=None, description="Minimum molecular weight (inclusive)"
    )
    mol_weight_max: float | None = strawberry.field(
        default=None, description="Maximum molecular weight (inclusive)"
    )
    ec_numbers: list[str] | None = None


@strawberry.input(description="Filter reactions by any field")
class ReactionFilter:
    """Filter input for querying reactions by any field.

    All fields are optional. Multiple fields are combined with AND.
    For range queries (e.g., seqLength), use min/max variants.
    """

    reversible: bool | None = None
