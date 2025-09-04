"""Data models for the Pyeed system using Pydantic v2."""

import re
import warnings
from dataclasses import dataclass
from enum import Enum
from typing import Annotated, Any, Dict, List, Optional
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field, PrivateAttr, field_validator


@dataclass(frozen=True)
class NodeHint:
    """Neo4j-specific metadata for model fields."""

    unique: bool = False
    index: bool = False
    vector_index: bool = False


@dataclass(frozen=True)
class EdgeHint:
    """Neo4j-specific metadata for model fields."""

    name: str
    outgoing_class_name: str  # Use string instead of actual class type
    outgoing_attr_name: str


class AnnotationType(str, Enum):
    """Protein and DNA annotation types."""

    ACTIVE_SITE = "active_site"
    ALLOSTERIC_SITE = "allosteric_site"
    ALPHAHELIX = "alpha_helix"
    BETASTRAND = "beta_strand"
    BINDING_SITE = "binding_site"
    MATURE_PROTEIN = "mature_protein"
    CODING_SEQ = "coding_sequence"
    DNA = "DNA"
    DOMAIN = "domain"
    FAMILY = "family"
    MOTIVE = "motive"
    PROTEIN = "protein"
    TURN = "turn"
    SIGNAL = "signal"
    PROPEP = "propeptide"


class BaseNode(BaseModel):
    """Base class for all nodes in the system."""

    model_config = ConfigDict(frozen=False, validate_assignment=True)

    custom: dict[str, Any] = Field(
        default_factory=dict, description="Arbitrary custom data as key-value pairs"
    )
    _id: str = PrivateAttr(default_factory=lambda: str(uuid4()))

    @field_validator("custom")
    @classmethod
    def validate_custom_keys(cls, v: dict[str, Any], info: Any) -> dict[str, Any]:
        """Validate that custom keys don't conflict with existing attributes."""
        if not v:
            return v

        # Get all field names from the current class and its parents
        field_names = set()
        current_class = info.data.get("__class__", cls)

        # Collect field names from current class and all parent classes
        while current_class and current_class != BaseModel:
            field_names.update(current_class.model_fields.keys())
            current_class = (
                current_class.__bases__[0] if current_class.__bases__ else None
            )

        # Check for conflicts
        conflicting_keys = []
        invalid_keys = []

        # Regex pattern for valid Python variable names
        valid_var_pattern = re.compile(r"^[a-zA-Z_][a-zA-Z0-9_]*$")

        for key, value in v.items():
            # Check for conflicts with existing attributes
            if key in field_names:
                conflicting_keys.append(key)
            elif key == "custom":
                conflicting_keys.append(key)

            # Check if key is a valid Python variable name
            if not valid_var_pattern.match(key):
                invalid_keys.append(key)

            # Check for nested dictionaries (not allowed)
            if isinstance(value, dict):
                raise ValueError(
                    f"Nested dictionaries are not allowed in custom fields. "
                    f"Key '{key}' contains a dictionary value."
                )

        if conflicting_keys:
            raise ValueError(
                f"Custom field keys cannot conflict with existing attributes or be 'custom': {conflicting_keys}"
            )

        if invalid_keys:
            raise ValueError(
                f"Invalid custom field keys: {invalid_keys}. "
                f"Keys must start with a letter or underscore and contain only letters, digits, or underscores."
            )

        return v

    def get_unique_id(self) -> str:
        """Return the value of the field marked with NodeHint(unique=True), else fallback to _id."""
        for field_name, field_info in type(self).model_fields.items():
            # Check for NodeHint metadata in Annotated fields
            for meta in getattr(field_info, "metadata", ()):
                if isinstance(meta, NodeHint) and meta.unique:
                    val = getattr(self, field_name, None)
                    if val is not None:
                        return str(val)
                    break

        return self._id

    def get_unique_field_name(self) -> str:
        """Return the name of the field marked with NodeHint(unique=True), else fallback to '_id'."""
        for field_name, field_info in type(self).model_fields.items():
            # Check for NodeHint metadata in Annotated fields
            for meta in getattr(field_info, "metadata", ()):
                if isinstance(meta, NodeHint) and meta.unique:
                    return field_name

        return "_id"

    @classmethod
    def get_unique_field_name_for_class(cls) -> str:
        """Return the name of the unique field for this class."""
        for field_name, field_info in cls.model_fields.items():
            # Check for NodeHint metadata in Annotated fields
            for meta in getattr(field_info, "metadata", ()):
                if isinstance(meta, NodeHint) and meta.unique:
                    return field_name

        return "_id"


class Organism(BaseNode):
    """Organism information."""

    accession_id: Annotated[
        str,
        EdgeHint(
            name="ORIGINATES_FROM",
            outgoing_class_name="Protein",
            outgoing_attr_name="accession_id",
        ),
    ] = Field(
        ...,
        description="Organism accession identifier",
    )
    tax_id: Annotated[int, NodeHint(unique=True)] = Field(
        ...,
        description="NCBI taxonomy ID",
    )
    name: Optional[str] = Field(
        None,
        description="Organism name",
    )

    @field_validator("tax_id")
    @classmethod
    def validate_tax_id(cls, v: int) -> int:
        if v <= 0:
            raise ValueError("Taxonomy ID must be positive")
        return v


class SequenceAnnotation(BaseNode):
    """Sequence annotation with positions and metadata."""

    accession_id: Annotated[
        str,
        EdgeHint(
            name="HAS_ANNOTATION",
            outgoing_class_name="Protein",
            outgoing_attr_name="accession_id",
        ),
    ] = Field(
        ...,
        description="Protein accession identifier",
    )
    annotation_type: AnnotationType = Field(
        ...,
        description="Type of annotation",
    )
    positions: List[int] = Field(
        ...,
        description="Sorted list of positions",
    )

    @field_validator("positions")
    @classmethod
    def validate_positions(cls, v: List[int]) -> List[int]:
        if not v:
            raise ValueError("Positions cannot be empty")
        return sorted(set(v))


class Molecule(BaseNode):
    """Chemical molecule information."""

    chebi_id: Annotated[str, NodeHint(unique=True)] = Field(
        ...,
        description="ChEBI identifier",
    )
    rhea_compound_id: Optional[str] = Field(
        None,
        description="RHEA compound identifier",
    )
    smiles: Optional[str] = Field(
        None,
        description="SMILES representation",
    )
    embedding: Optional[List[float]] = Field(
        None,
        description="Embedding vector",
    )


class Reaction(BaseNode):
    """Chemical reaction information."""

    accession_id: Annotated[
        str,
        EdgeHint(
            name="HAS_REACTION",
            outgoing_class_name="Protein",
            outgoing_attr_name="accession_id",
        ),
    ] = Field(
        ...,
        description="Protein accession identifier that can catalyze the reaction",
    )
    rhea_id: str = Field(
        ...,
        description="RHEA reaction identifier",
    )
    substrate: Annotated[
        List[Molecule],
        EdgeHint(
            name="HAS_SUBSTRATE",
            outgoing_class_name="Molecule",
            outgoing_attr_name="accession_id",
        ),
    ] = Field(
        default_factory=list,
        description="List of ChEBI identifiers",
    )
    products: List[Molecule] = Field(
        default_factory=list,
        description="List of ChEBI identifiers",
    )
    reversible: bool = Field(
        default=False,
        description="Whether the reaction is reversible",
    )


class GOAnnotation(BaseNode):
    """Gene Ontology annotation."""

    go_id: Annotated[str, NodeHint(unique=True)] = Field(
        ...,
        description="Gene Ontology identifier",
    )
    term: Optional[str] = Field(
        None,
        description="GO term name",
    )
    definition: Optional[str] = Field(
        None,
        description="GO term definition",
    )

    # For GO annotations, we need a way to link them to proteins
    # This will be used when creating GO-Protein relationships
    protein_accession_id: Annotated[
        Optional[str],
        EdgeHint(
            name="HAS_GO_ANNOTATION",
            outgoing_class_name="Protein",
            outgoing_attr_name="accession_id",
        ),
    ] = Field(
        None,
        description="Protein accession identifier (for relationship creation)",
    )


class Embedding(BaseNode):
    """Metadata about an embedding."""

    accession_id: Annotated[
        str,
        EdgeHint(
            name="HAS_EMBEDDING",
            outgoing_class_name="Protein",
            outgoing_attr_name="accession_id",
        ),
    ] = Field(
        ...,
        description="Protein accession identifier",
    )
    description: str = Field(
        ...,
        description="Description of the embedding",
    )
    model_name: str = Field(
        ...,
        description="Name of the embedding model",
    )
    layer_index: int = Field(
        default=-1,
        description="Model layer number embedding matrix was extracted prior to pooling.",
    )
    pooling_method: str = Field(
        ...,
        description="Pooling method",
    )
    vector: Annotated[
        List[float],
        NodeHint(vector_index=True),
    ] = Field(
        ...,
        description="The embedding vector",
    )
    n_dims: int = Field(
        ...,
        description="Embedding vector length",
    )


class Protein(BaseNode):
    """Protein sequence and metadata."""

    accession_id: Annotated[
        str,
        NodeHint(unique=True),
    ] = Field(
        ...,
        description="Protein accession identifier",
    )
    sequence: str = Field(
        ...,
        description="Amino acid sequence",
    )
    name: Optional[str] = Field(
        None,
        description="Protein name",
    )
    seq_length: Optional[int] = Field(
        None,
        description="Sequence length",
    )
    organisms: List[Organism] = Field(
        default_factory=list,
        description="Organisms the protein originates from",
    )
    mol_weight: Optional[float] = Field(
        None,
        description="Molecular weight in Daltons",
    )
    ec_number: Optional[str] = Field(
        None,
        description="Enzyme Commission number",
        pattern="^\d{1,5}\.\d{1,5}\.\d{1,5}\.\d{1,5}$",
    )
    nucleotide_id: Optional[str] = Field(
        None,
        description="Associated nucleotide ID",
    )
    nucleotide_start: Optional[int] = Field(
        None,
        description="Nucleotide start position",
    )
    nucleotide_end: Optional[int] = Field(
        None,
        description="Nucleotide end position",
    )
    locus_tag: Optional[str] = Field(
        None,
        description="Locus tag",
    )
    structure_ids: List[str] = Field(
        default_factory=list,
        description="Structure identifiers",
    )
    go_terms: List[GOAnnotation] = Field(
        default_factory=list,
        description="GO term identifiers",
    )
    reactions: List[Reaction] = Field(
        default_factory=list,
        description="RHEA reaction identifiers",
    )
    embeddings: Dict[str, Embedding] = Field(
        default_factory=dict,
        description="Embeddings from different models",
    )
    annotations: Dict[str, SequenceAnnotation] = Field(
        default_factory=dict,
        description="Sequence annotations",
    )

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("Sequence cannot be empty")
        if not all(c in "ACDEFGHIKLMNPQRSTVWY" for c in v.upper()):
            warnings.warn(
                f"Sequence contains non-standard characters: {set(v) - set('ACDEFGHIKLMNPQRSTVWY')}"
            )
        return v.upper()

    @field_validator("seq_length")
    @classmethod
    def validate_seq_length(cls, v: Optional[int], info: Any) -> Optional[int]:
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
                f"Provided sequence length {v} does not match actual sequence length {actual_length}"
            )
        elif v <= 0:
            raise ValueError("Sequence length must be positive")

        return v

    @field_validator("annotations")
    @classmethod
    def validate_annotations(
        cls, v: Dict[str, SequenceAnnotation]
    ) -> Dict[str, SequenceAnnotation]:
        """Ensure each annotation type occurs only once."""
        seen_types = set()
        for annotation in v.values():
            if annotation.annotation_type in seen_types:
                raise ValueError(
                    f"Annotation type {annotation.annotation_type} appears multiple times"
                )
            seen_types.add(annotation.annotation_type)
        return v

    def add_embedding(self, embedding: Embedding, replace: bool = False) -> None:
        """Add or replace embedding for a specific model.

        Args:
            embedding: Embedding object to add.
            replace: If False, raises if embedding exists; if True, replaces existing.

        Raises:
            ValueError: If embedding exists and replace is False.
        """
        if embedding.accession_id != self.accession_id:
            raise ValueError(
                f"Embedding accession_id {embedding.accession_id} does not match protein accession_id {self.accession_id}"
            )

        if embedding.description in self.embeddings and not replace:
            raise ValueError(
                f"Embedding for model '{embedding.description}' already exists. "
                "Set replace=True to overwrite."
            )
        self.embeddings[embedding.description] = embedding

    def remove_embedding(self, model_name: str) -> None:
        """Remove embedding for a specific model."""
        self.embeddings.pop(model_name, None)

    def add_annotation(
        self, annotation: SequenceAnnotation, replace: bool = False
    ) -> None:
        """Add or replace annotation for a specific type.

        Args:
                    annotation: SequenceAnnotation object to add.
                    replace: If False, raises if annotation exists; if True, replaces existing.

                Raises:
                    ValueError: If annotation exists and replace is False.
        """
        # check that accession_id is the same as the annotation
        if annotation.accession_id != self.accession_id:
            raise ValueError(
                f"Annotation accession_id {annotation.accession_id} does not match protein accession_id {self.accession_id}"
            )

        annotation_type = annotation.annotation_type.value

        if annotation_type in self.annotations and not replace:
            raise ValueError(
                f"Annotation for type '{annotation_type}' already exists. "
                "Set replace=True to overwrite."
            )

        self.annotations[annotation_type] = annotation

    def remove_annotation(self, annotation_type: AnnotationType) -> None:
        """Remove annotation of specified type."""
        self.annotations.pop(annotation_type.value, None)

    def get_annotation(
        self, annotation_type: AnnotationType
    ) -> Optional[SequenceAnnotation]:
        """Get annotation of specified type."""
        return self.annotations.get(annotation_type.value)

    def has_annotation_type(self, annotation_type: AnnotationType) -> bool:
        """Check if protein has annotation of specified type."""
        return annotation_type.value in self.annotations

    def add_organism(self, organism: Organism) -> None:
        """Add organism to protein."""
        if organism.accession_id != self.accession_id:
            raise ValueError(
                f"Organism accession_id {organism.accession_id} does not match protein accession_id {self.accession_id}"
            )
        if organism.tax_id in [o.tax_id for o in self.organisms]:
            raise ValueError(
                f"Organism with tax_id {organism.tax_id} already exists in protein {self.accession_id}"
            )
        self.organisms.append(organism)

    def remove_organism(self, organism: Organism) -> None:
        """Remove organism from protein."""
        self.organisms.remove(organism)


class DNA(BaseNode):
    """DNA sequence and metadata."""

    accession_id: str = Field(
        ...,
        description="DNA accession identifier",
    )
    sequence: str = Field(
        ...,
        description="DNA sequence",
    )
    name: Optional[str] = Field(
        None,
        description="DNA name",
    )
    seq_length: int = Field(
        ...,
        description="Sequence length",
    )
    go_terms: List[str] = Field(
        default_factory=list,
        description="GO term identifiers",
    )
    embedding: List[float] = Field(
        default_factory=list,
        description="DNA embedding vector",
    )
    gc_content: Optional[float] = Field(
        None,
        description="GC content percentage",
    )

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("Sequence cannot be empty")
        if not all(c in "ACGTN" for c in v):
            warnings.warn(
                f"Sequence contains non-standard characters: {set(v) - set('ACGTN')}"
            )
        return v

    @field_validator("gc_content")
    @classmethod
    def validate_gc_content(cls, v: Optional[float]) -> Optional[float]:
        if v is not None and (v < 0 or v > 100):
            raise ValueError("GC content must be between 0 and 100")
        return v


# Example usage
if __name__ == "__main__":
    from rich import print

    prot = Protein(
        accession_id="P01234",
        sequence="MALWMRLLPLLALLALWGPDPAAA",
        name="Test Protein",
        seq_length=24,
        mol_weight=1000,
        ec_number="1.2.3.4",
        nucleotide_id="N01234",
        nucleotide_start=1,
        nucleotide_end=25,
        locus_tag="L01234",
        custom={
            "mentioned_in": [
                "doi:10.1016/j.xinn.2025.100344",
                "doi:10.1016/j.xinn.2025.100345",
            ],
            "already_characterized": False,
            "dictd": 123,
        },
    )

    prot.add_annotation(
        SequenceAnnotation(
            accession_id="P01234",
            annotation_type=AnnotationType.ACTIVE_SITE,
            positions=[1, 2, 3],
            custom={"my_custom_evidence": "literature", "validated": True},
        )
    )

    # add second annotation
    prot.add_annotation(
        SequenceAnnotation(
            accession_id="P01234",
            annotation_type=AnnotationType.BINDING_SITE,
            positions=[4, 5, 6],
        )
    )

    # add embedding
    prot.add_embedding(
        Embedding(
            accession_id="P01234",
            description="esm2_mean_pooled",
            model_name="esm2-t33-650M-UR50S",
            pooling_method="mean",
            vector=[0.1, 0.2, 0.3],
            n_dims=3,
        ),
    )

    # add another embedding
    prot.add_embedding(
        Embedding(
            accession_id="P01234",
            description="esm2_mean_pooled_last_hidden_state",
            model_name="esm2-t33-650M-UR50S",
            pooling_method="mean",
            vector=[0.4, 0.5, 0.6],
            n_dims=3,
        ),
    )
    prot.add_embedding(
        Embedding(
            accession_id="P01234",
            description="esm2_range_pooled",
            model_name="esm2-t33-650M-UR50S",
            pooling_method="range",
            vector=[0.4, 0.5, 0.6, 0.7, 0.8],
            n_dims=5,
        ),
    )

    print(prot)

    for embedding in prot.embeddings.values():
        print(embedding._id)
