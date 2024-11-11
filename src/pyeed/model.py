from enum import Enum

from neomodel import (
    ArrayProperty,
    FloatProperty,
    IntegerProperty,
    RelationshipTo,
    StringProperty,
    StructuredRel,
    UniqueIdProperty,
    VectorIndex,
    UniqueProperty,
    StructuredNode,
)

# from pyeed.nodes_and_relations import StrictStructuredNode

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
    MATURE_PROTEIN = "mature protein"
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


class StandardNumberingRel(StructuredRel):
    positions = ArrayProperty(IntegerProperty(), required=True)

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: StrictStructuredNode,
        molecule2: StrictStructuredNode,
        positions: list[str],
    ):
        molecule1.sequences_protein.connect(
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


class StandardNumbering(StrictStructuredNode):
    name = StringProperty(required=True)
    definition = StringProperty(required=True)

    # Relationships
    sequences_protein = RelationshipTo(
        "Protein", "HAS_STANDARD_NUMBERING", model=StandardNumberingRel
    )


class PairwiseAlignmentResult(StructuredRel):
    """A relationship representing the similarity between two sequences.

    Args:
        similarity (float): The similarity score between the two sequences.
        gaps (int): The number of gaps in the alignment.
        mismatches (int): The number of mismatches in the alignment.
        score (int): The alignment score.
        query_aligned (str): The aligned sequence of the query.
        target_aligned (str): The aligned sequence of the target.
    """

    similarity = FloatProperty(required=True)
    gaps = IntegerProperty(required=True)
    mismatches = IntegerProperty(required=True)
    score = IntegerProperty()
    query_aligned = StringProperty()
    target_aligned = StringProperty()

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: StrictStructuredNode,
        molecule2: StrictStructuredNode,
        similarity: float,
        gaps: int,
        mismatches: int,
        score: int,
        query_aligned: str,
        target_aligned: str,
    ):
        """Validates the similarity and connects the two molecules.

        Args:
            molecule1 (StrictStructuredNode): Protein or DNA node
            molecule2 (StrictStructuredNode): Protein or DNA node
            similarity (float): Percentage similarity between the two sequences.
            gaps (int): Number of gaps in the alignment.
            mismatches (int): Number of mismatches in the alignment.
            score (int): Alignment score.
            query_aligned (str): Alignment of the query sequence.
            target_aligned (str): Alignment of the target sequence.

        Returns:
            Similarity: The created similarity relationship.
        """
        molecule1.similar.connect(
            molecule2,
            {
                "similarity": similarity,
                "gaps": gaps,
                "mismatches": mismatches,
                "score": score,
                "query_aligned": query_aligned,
                "target_aligned": target_aligned,
            },
        )

        return cls(
            similarity=similarity,
            gaps=gaps,
            mismatches=mismatches,
            score=score,
            query_aligned=query_aligned,
            target_aligned=target_aligned,
        )


class GOAnnotation(StrictStructuredNode):
    """A Gene Ontology annotation for a protein or dna sequence.

    Args:
        go_id (str): The Gene Ontology ID.
        term (str): The name of the GO term.
        definition (str): The definition of the GO term.
    """

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
        from_monomer (str): The original monomer at the mutation position.
        to_monomer (str): The mutated residue at the mutation
    """

    from_position = IntegerProperty(required=True)
    to_position = IntegerProperty(required=True)
    from_monomer = StringProperty(required=True)
    to_monomer = StringProperty(required=True)

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

        Returns:
            Mutation: The created mutation relationship.

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
                "from_monomer": from_monomer,
                "to_monomer": to_monomer,
            },
        )

        return cls(
            from_position=from_position,
            to_position=to_position,
            from_monomer=from_monomer,
            to_monomer=to_monomer,
        )

    @property
    def label(self):
        return f"{self.from_monomer}{self.from_position}{self.to_monomer}"


class Protein(StrictStructuredNode):
    """A protein sequence node in the database."""

    accession_id = StringProperty(unique_index=True, required=True)
    sequence = StringProperty(required=True)
    name = StringProperty()
    seq_length = IntegerProperty(required=True)
    mol_weight = FloatProperty()
    ec_number = StringProperty()
    nucleotide_id = StringProperty()
    nucleotide_start = IntegerProperty()
    nucleotide_end = IntegerProperty()
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
    pairwise_aligned = RelationshipTo(
        "Protein", "PAIRWISE_ALIGNED", model=PairwiseAlignmentResult
    )


class DNA(StrictStructuredNode):
    """A DNA sequence node in the database."""

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
    pairwise_aligned = RelationshipTo(
        "DNA", "PAIRWISE_ALIGNED", model=PairwiseAlignmentResult
    )
