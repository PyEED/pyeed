import re
import warnings
from dataclasses import dataclass
from enum import Enum
from typing import Annotated, Any, Dict, Iterable, List, Optional, Tuple, Type
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field, ValidationInfo, field_validator
from typing_extensions import get_args, get_origin


@dataclass(frozen=True)
class NodeHint:
    unique: bool = False
    index: bool = False
    vector_index: bool = False


def node_label(cls: Type[BaseModel]) -> str:
    return cls.__name__


def extract_label_key_map(models: Iterable[Type[BaseModel]]) -> dict[str, str]:
    """{Label: unique_key} (first field marked unique=True per model)."""
    out: dict[str, str] = {}
    for m in models:
        for name, ann in m.model_fields.items():  # pydantic v2
            tp = ann.annotation
            if get_origin(tp) is Annotated:
                base, *meta = get_args(tp)
                for x in meta:
                    if isinstance(x, NodeHint) and x.unique:
                        out[node_label(m)] = name
                        break
                if node_label(m) in out:
                    break
    return out


def extract_indexes(
    models: Iterable[Type[BaseModel]],
) -> list[tuple[str, str, NodeHint]]:
    """[(Label, prop, hint)] for index/vector_index flags."""
    out = []
    for m in models:
        lbl = node_label(m)
        for name, ann in m.model_fields.items():
            tp = ann.annotation
            if get_origin(tp) is Annotated:
                base, *meta = get_args(tp)
                for x in meta:
                    if isinstance(x, NodeHint) and (x.index or x.vector_index):
                        out.append((lbl, name, x))
    return out


class AnnotationType(str, Enum):
    """Protein and DNA annotation types."""

    SITE = "site"
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

    model_config = ConfigDict(
        frozen=False, validate_assignment=True, use_enum_values=True
    )

    custom: dict[str, Any] = Field(
        default_factory=dict, description="Arbitrary custom data as key-value pairs"
    )

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

    def to_dict(self) -> dict[str, Any]:
        """Flattens "custom" field to the top level of the dictionary"""
        nested_dict = self.model_dump(exclude_none=True, exclude_unset=True)
        custom = nested_dict.get("custom", {})
        if custom:
            nested_dict.pop("custom")
            return {**nested_dict, **custom}
        return nested_dict

    def get_unique_model_field(self) -> str:
        """Returns the name of the field marked with NodeHint(unique=True)"""
        print("getting unique model field for", type(self))
        for field_name, field_info in type(self).model_fields.items():
            print(field_name, field_info)
            if not field_info.metadata:
                continue
            if (
                isinstance(field_info.metadata[0], NodeHint)
                and field_info.metadata[0].unique
            ):
                print("found unique field", field_name)
                return field_name

        raise ValueError(
            f"No unique field found. No field of {type(self)} is marked with NodeHint(unique=True)"
        )

    def graphify(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Flatten this node and all nested BaseNode values into nodes & edges."""
        nodes: List[Dict[str, Any]] = []
        edges: List[Dict[str, Any]] = []

        def emit_node(obj: BaseNode) -> Dict[str, Any]:
            label = obj.__class__.__name__
            unique_field = obj.get_unique_model_field()
            unique_value = getattr(obj, unique_field)
            props = obj.to_dict()
            return {"label": label, "key": (unique_field, unique_value), "props": props}

        def recurse(parent: BaseNode, fname: str, val: Any) -> None:
            if isinstance(val, BaseNode):
                nodes.append(emit_node(val))
                edges.append(
                    {
                        "type": fname.upper(),
                        "src": (
                            parent.__class__.__name__,
                            parent.get_unique_model_field(),
                            getattr(parent, parent.get_unique_model_field()),
                        ),
                        "dst": (
                            val.__class__.__name__,
                            val.get_unique_model_field(),
                            getattr(val, val.get_unique_model_field()),
                        ),
                    }
                )
            elif isinstance(val, Iterable) and not isinstance(val, (str, bytes)):
                for item in val:
                    recurse(parent, fname, item)

        # Root node
        nodes.append(emit_node(self))

        # Process all fields
        for fname in self.model_dump().keys():
            recurse(self, fname, getattr(self, fname, None))

        return nodes, edges


class Organism(BaseNode):
    """Organism information."""

    tax_id: Annotated[int, NodeHint(unique=True)] = Field(
        ...,
        description="NCBI taxonomy ID",
    )
    name: Optional[str] = Field(
        default=None,
        description="Organism name",
    )

    @field_validator("tax_id")
    @classmethod
    def validate_tax_id(cls, v: int) -> int:
        if v <= 0:
            raise ValueError("Taxonomy ID must be positive")
        return v


class Annotation(BaseNode):
    """Sequence annotation with positions and metadata."""

    id: Annotated[str, NodeHint(unique=True)] = Field(
        default_factory=lambda: str(uuid4()),
        description="Annotation identifier",
    )
    annotation_type: AnnotationType = Field(
        ...,
        description="Type of annotation",
    )
    positions: List[int] = Field(
        ...,
        description="Sorted list of positions",
    )
    description: Optional[str] = Field(
        default=None,
        description="Description of the annotation",
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
    name: Optional[str] = Field(
        None,
        description="Molecule name",
    )
    smiles: Optional[str] = Field(
        None,
        description="SMILES representation",
    )
    inchi: Optional[str] = Field(
        None,
        description="InChI representation",
    )
    embedding: Optional[List[float]] = Field(
        None,
        description="Embedding vector",
    )


class Reaction(BaseNode):
    """Chemical reaction information."""

    rhea_id: Annotated[str, NodeHint(unique=True)] = Field(
        ...,
        description="RHEA reaction identifier",
    )
    description: Optional[str] = Field(
        None,
        description="Reaction description",
    )
    substrates: List[Molecule] = Field(
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


class Embedding(BaseNode):
    """Metadata about an embedding."""

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
    id: Annotated[str, NodeHint(unique=True)] = Field(
        default_factory=lambda: str(uuid4()),
        description="Embedding identifier",
    )

    @field_validator("model_name")
    @classmethod
    def model_name_slug(cls, v: str) -> str:
        # Convert to lowercase and replace non [a-z0-9_] with underscores
        v_clean = re.sub(r"[^a-z0-9_]", "_", v.lower())
        return v_clean

    @field_validator("pooling_method")
    @classmethod
    def pooling_method_slug(cls, v: str) -> str:
        # Convert to lowercase and replace non [a-z0-9_] with underscores
        v_clean = re.sub(r"[^a-z0-9_]", "_", v.lower())
        return v_clean

    @field_validator("vector")
    @classmethod
    def _check_vector_len(cls, v: List[float], info: ValidationInfo) -> List[float]:
        # Access n_dims via info.data (other fields that have already been validated)
        n_dims = info.data.get("n_dims")
        if n_dims is not None and n_dims != len(v):
            raise ValueError(f"n_dims={n_dims} != len(vector)={len(v)}")
        return v

    @property
    def neo4j_vector_prop(self) -> str:
        """Dynamic property to write into Neo4j (one ANN index per property)"""
        return f"vec__{self.model_name}__{self.pooling_method}"

    def to_dict(self) -> dict[str, Any]:
        """
        For DB write: flatten 'custom' onto root, and rename 'vector' to the dynamic
        property (vec__{model}__{pool}). Keep other fields as-is.
        """
        d = self.model_dump()
        custom = d.pop("custom", {}) or {}
        vector_value = d.pop("vector")
        vector_dict = {self.neo4j_vector_prop: vector_value}

        return {**d, **custom, **vector_dict}


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
    ec_numbers: Optional[List[str]] = Field(
        None,
        description="Enzyme Commission numbers associated with the protein",
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
    embeddings: List[Embedding] = Field(
        default_factory=list,
        description="Embeddings from different models",
    )
    annotations: List[Annotation] = Field(
        default_factory=list,
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
    def validate_annotations(cls, v: List[Annotation]) -> List[Annotation]:
        """Validate annotations list."""
        # No need for uniqueness validation - multiple annotations of same type are allowed
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
        ec_numbers=["1.2.3.4"],
        custom={
            "mentioned_in": [
                "doi:10.1016/j.xinn.2025.100344",
                "doi:10.1016/j.xinn.2025.100345",
            ],
            "already_characterized": False,
            "dictd": 123,
        },
    )

    prot.annotations.append(
        Annotation(
            annotation_type=AnnotationType.ACTIVE_SITE,
            positions=[1, 2, 3],
            custom={"my_custom_evidence": "literature", "validated": True},
        ),
    )

    # add second annotation
    prot.annotations.append(
        Annotation(
            annotation_type=AnnotationType.BINDING_SITE,
            positions=[4, 5, 6],
        ),
    )

    # add embedding
    prot.embeddings.append(
        Embedding(
            model_name="esm2_t33_650M_UR50S",
            pooling_method="mean",
            vector=[0.1, 0.2, 0.3],
            n_dims=3,
        ),
    )

    # add another embedding
    prot.embeddings.append(
        Embedding(
            model_name="esm2_t33_650M_UR50S",
            pooling_method="mean",
            vector=[0.4, 0.5, 0.6],
            n_dims=3,
        ),
    )
    prot.embeddings.append(
        Embedding(
            model_name="esm2-t33-650M-UR50S",
            pooling_method="range",
            vector=[0.4, 0.5, 0.6, 0.7, 0.8],
            n_dims=5,
        ),
    )

    print(prot.graphify())
