from enum import Enum

from neomodel import (
    ArrayProperty,
    FloatProperty,
    IntegerProperty,
    RelationshipTo,
    StringProperty,
    StructuredNode,
    StructuredRel,
    UniqueIdProperty,
    UniqueProperty,
    VectorIndex,
)


class StrictStructuredNode(StructuredNode):
    """A StructuredNode subclass that raises an error if an invalid property is provided."""

    __abstract_node__ = True

    def __init__(self, *args, **kwargs):
        # Get the defined properties of the model
        allowed_properties = set(self.__class__._class_properties())

        # Check if any provided properties are not in the allowed set
        for key in kwargs:
            if key not in allowed_properties:
                raise AttributeError(
                    f"'{key}' is not a valid property for {self.__class__.__name__}"
                )

        super().__init__(*args, **kwargs)

    @classmethod
    def _class_properties(cls):
        """Retrieve all allowed properties (fields) defined on the class."""
        return {
            k
            for k, v in cls.__dict__.items()
            if isinstance(
                v,
                (
                    StringProperty,
                    IntegerProperty,
                    FloatProperty,
                    ArrayProperty,
                    UniqueIdProperty,
                ),
            )
        }

    def save(self, *args, **kwargs):
        """Validates the properties and then saves the node."""
        allowed_properties = self.__class__._class_properties()

        # Only validate properties defined in the model schema
        for field, prop in self.__dict__.items():
            if field not in allowed_properties:
                continue  # Skip non-class properties (like internal Neo4j fields)

            if prop is None or callable(prop):
                continue

            try:
                neo_type = getattr(self.__class__, field)
            except AttributeError:
                raise AttributeError(
                    f"'{self.__class__.__name__}' has no attribute '{field}'"
                )

            # Skip validation for UniqueIdProperty
            if isinstance(neo_type, UniqueIdProperty):
                continue

            # Validate StringProperty
            if isinstance(neo_type, StringProperty) and not isinstance(prop, str):
                raise TypeError(
                    f"Expected a string for '{field}', got {type(prop).__name__}"
                )

            # Validate IntegerProperty
            elif isinstance(neo_type, IntegerProperty) and not isinstance(prop, int):
                raise TypeError(
                    f"Expected an integer for '{field}', got {type(prop).__name__}"
                )

            # Validate FloatProperty
            elif isinstance(neo_type, FloatProperty) and not isinstance(prop, float):
                raise TypeError(
                    f"Expected a float for '{field}', got {type(prop).__name__}"
                )

            # Validate ArrayProperty
            elif isinstance(neo_type, ArrayProperty):
                if not isinstance(prop, list):
                    raise TypeError(
                        f"Expected a list for '{field}', got {type(prop).__name__}"
                    )

                # Validate list of integers, strings, or floats
                base_property = neo_type.base_property
                if isinstance(base_property, StringProperty):
                    if not all(isinstance(item, str) for item in prop):
                        raise TypeError(f"All items in '{field}' must be strings")
                elif isinstance(base_property, IntegerProperty):
                    if not all(isinstance(item, int) for item in prop):
                        raise TypeError(f"All items in '{field}' must be integers")
                elif isinstance(base_property, FloatProperty):
                    if not all(isinstance(item, float) for item in prop):
                        raise TypeError(f"All items in '{field}' must be floats")

        return super().save(*args, **kwargs)

    @classmethod
    def get_or_save(cls, **kwargs):
        """Attempts to save the node first, and if it already exists (due to unique constraint), retrieves it."""
        try:
            # Attempt to create and save a new node
            instance = cls(**kwargs)
            instance.save()
            return instance
        except UniqueProperty:
            # If a unique constraint error occurs, retrieve the existing node
            return cls.nodes.get(**kwargs)



class Annotation(Enum):
    ACTIVE_SITE = "active site"
    ALLOSTERIC_SITE = "allosteric site"
    ALPHAHELIX = "alpha helix"
    BETASTRAND = "beta strand"
    BINDING_SITE = "binding site"
    CODING_SEQ = "coding sequence"
    DNA = "DNA"
    DOMAIN = "domain"
    FAMILY = "family"
    MOTIVE = "motive"
    PROTEIN = "protein"

class Organism(StrictStructuredNode):
    taxonomy_id: int = IntegerProperty(required=True, unique_index=True)
    name = StringProperty()

class SiteRel(StructuredRel):
    positions = ArrayProperty(IntegerProperty(), required=True)

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: StrictStructuredNode,
        molecule2: StrictStructuredNode,
        positions: list,
    ):
        molecule1.site.connect(
            molecule2,
            {
                "positions": positions,
            },
        )

        return cls(
            positions=positions,
        )

    @property
    def label(self):
        return f"{self.positions}"

class Site(StrictStructuredNode):
    site_id = UniqueIdProperty()
    name = StringProperty()
    annotation = StringProperty(
        choices=[(e.value, e.name) for e in Annotation], required=True
    )

class Region(StrictStructuredNode):
    region_id = UniqueIdProperty()
    # start = IntegerProperty(required=True)
    # end = IntegerProperty(required=True)
    annotation = StringProperty(
        choices=[(e.value, e.name) for e in Annotation], required=True
    )

class RegionRel(StructuredRel):
    start = IntegerProperty(required=True)
    end = IntegerProperty(required=True)

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: StrictStructuredNode,
        molecule2: StrictStructuredNode,
        start: int,
        end: int,
    ):
        molecule1.region.connect(
            molecule2,
            {
                "start": start,
                "end": end,
            },
        )

        return cls(
            start=start,
            end=end,
        )

    @property
    def label(self):
        return f"{self.start}-{self.end}"

class Similarity(StructuredRel):
    similarity = FloatProperty(required=True)
    gaps = IntegerProperty()
    mismatches = IntegerProperty()
    score = IntegerProperty()
    aligned_sequences = ArrayProperty(StringProperty())

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: StrictStructuredNode,
        molecule2: StrictStructuredNode,
        similarity: float,
        gaps: int,
        mismatches: int,
        score: int,
        aligned_sequences: list,
    ):
        molecule1.similar.connect(
            molecule2,
            {
                "similarity": similarity,
                "gaps": gaps,
                "mismatches": mismatches,
                "score": score,
                "aligned_sequences": aligned_sequences,
            },
        )

        return cls(
            similarity=similarity,
            gaps=gaps,
            mismatches=mismatches,
            score=score,
            aligned_sequences=aligned_sequences,
        )

class GOAnnotation(StrictStructuredNode):
    go_id = StringProperty(unique_index=True, required=True)
    term = StringProperty()
    definition = StringProperty()

    @property
    def label(self):
        return self.term


class Mutation(StructuredRel):
    """A relationship representing a mutation between two sequences.

    Args:
        from_position (int): The position of the mutation in the original sequence.
        to_position (int): The position of the mutation in the mutated sequence.
        from_residue (str): The original residue at the mutation position.
        to_residue (str): The mutated residue at the mutation
    """

    from_position = IntegerProperty(required=True)
    to_position = IntegerProperty(required=True)
    from_residue = StringProperty(required=True)
    to_residue = StringProperty(required=True)

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: StrictStructuredNode,
        molecule2: StrictStructuredNode,
        from_position: int,
        to_position: int,
        from_monomer: str,
        to_monomer: str,
    ):
        """Validates the mutation and connects the two molecules.

        Args:
            molecule1 (StrictStructuredNode): DNA or Protein node
            molecule2 (StrictStructuredNode): DNA or Protein node
            from_position (int): Position of the mutation in the original sequence. 0-indexed.
            to_position (int): Position of the mutation in the mutated sequence. 0-indexed.
            from_monomer (str): Original residue / nucleotide at the specified position.
            to_monomer (str): Mutated residue / nucleotide at the specified position.

        Raises:
            ValueError: If the specified positions or residues do not match the sequences.
        """

        if molecule1.sequence[from_position] != from_monomer:
            raise ValueError(
                f"Monomer '{from_monomer}' does not match the sequence {molecule1.accession_id} at position {from_position}"
            )

        if molecule2.sequence[to_position] != to_monomer:
            raise ValueError(
                f"Monomer '{to_monomer}' does not match the sequence {molecule2.accession_id} at position {to_position}"
            )

        molecule1.mutation.connect(
            molecule2,
            {
                "from_position": from_position,
                "to_position": to_position,
                "from_residue": from_monomer,
                "to_residue": to_monomer,
            },
        )

        return cls(
            from_position=from_position,
            to_position=to_position,
            from_residue=from_monomer,
            to_residue=to_monomer,
        )

    @property
    def label(self):
        return f"{self.from_residue}{self.from_position}{self.to_residue}"


class Protein(StrictStructuredNode):
    accession_id = StringProperty(unique_index=True, required=True)
    sequence = StringProperty(required=True)
    name = StringProperty()
    seq_length = IntegerProperty(required=True)
    mol_weight = FloatProperty()
    ec_number = StringProperty()
    nucleotide_id = StringProperty()
    locus_tag = StringProperty()
    structure_ids = ArrayProperty(StringProperty())
    go_terms = ArrayProperty(StringProperty())
    embedding = ArrayProperty(
        FloatProperty(),
        vector_index=VectorIndex(dimensions=1280),
    )

    # Relationships
    organism = RelationshipTo("Organism", "ORIGINATES_FROM")
    site = RelationshipTo("Site", "HAS_SITE", model=SiteRel)
    region = RelationshipTo("Region", "HAS_REGION", model=RegionRel)
    go_annotation = RelationshipTo("GOAnnotation", "ASSOCIATED_WITH")
    mutation = RelationshipTo("Protein", "MUTATION", model=Mutation)
    coding_sequence = RelationshipTo("DNA", "HAS_CODING_SEQUENCE")
    similar = RelationshipTo("Protein", "SIMILARITY", model=Similarity)


class DNA(StrictStructuredNode):
    accession_id = StringProperty(unique_index=True, required=True)
    sequence = StringProperty(required=True)
    name = StringProperty()
    seq_length = IntegerProperty(required=True)
    go_terms = ArrayProperty(StringProperty())
    embedding = ArrayProperty(
        FloatProperty(),
        vector_index=VectorIndex(dimensions=1048),
    )
    gc_content = FloatProperty()

    # Relationships
    organism = RelationshipTo("Organism", "ORIGINATES_FROM")
    site = RelationshipTo("Site", "HAS_SITE", model=SiteRel)
    region = RelationshipTo("Region", "HAS_REGION", model=RegionRel)
    go_annotation = RelationshipTo("GOAnnotation", "ASSOCIATED_WITH")
    mutation = RelationshipTo("DNA", "MUTATION", model=Mutation)
    protein = RelationshipTo("Protein", "ENCODES", model=RegionRel)
    similar = RelationshipTo("DNA", "SIMILARITY", model=Similarity)
