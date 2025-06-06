from enum import Enum
from typing import Any, cast

# from pyeed.nodes_and_relations import StrictStructuredNode
from neomodel import (
    ArrayProperty,
    BooleanProperty,
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


class StrictStructuredNode(StructuredNode):  # type: ignore
    """A StructuredNode subclass that raises an error if an invalid property is provided."""

    __abstract_node__ = True

    def __init__(self, *args: Any, **kwargs: Any) -> None:
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
    def _class_properties(cls) -> set[str]:
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

    def save(self, *args: Any, **kwargs: Any) -> None:
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

            # Validate BoleanProperty
            elif isinstance(neo_type, BooleanProperty) and not isinstance(prop, bool):
                raise TypeError(
                    f"Expected a boolean for '{field}', got {type(prop).__name__}"
                )

        super().save(*args, **kwargs)  # Don't return the result

    @classmethod
    def get_or_save(cls, **kwargs: Any) -> StructuredNode:
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
    TURN = "turn"
    SIGNAL = "signal"
    PROPEP = "propep"


class Organism(StrictStructuredNode):
    taxonomy_id = IntegerProperty(required=True, unique_index=True)
    name = StringProperty()

    @classmethod
    def get_or_save(cls, **kwargs: Any) -> "Organism":
        taxonomy_id = kwargs.get("taxonomy_id")
        name = kwargs.get("name")
        try:
            organism = cast(Organism, cls.nodes.get(taxonomy_id=taxonomy_id))
            return organism
        except cls.DoesNotExist:
            try:
                organism = cls(taxonomy_id=taxonomy_id, name=name)
                organism.save()
                return organism
            except Exception as e:
                print(f"Error during saving of the organism: {e}")
                raise


class Mutation(StructuredRel):  # type: ignore
    """A relationship representing mutations between two sequences."""

    from_positions = ArrayProperty(IntegerProperty(), required=True)
    to_positions = ArrayProperty(IntegerProperty(), required=True)
    from_monomers = ArrayProperty(StringProperty(), required=True)
    to_monomers = ArrayProperty(StringProperty(), required=True)

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: Any,
        molecule2: Any,
        from_positions: list[int],
        to_positions: list[int],
        from_monomers: list[str],
        to_monomers: list[str],
    ) -> "Mutation":
        """Validates the mutations and connects the two molecules, ensuring that no double mutations
        occur – i.e. if a mutation affecting any of the same positions already exists between these proteins,
        a new mutation cannot be created.

        Raises:
            ValueError: If input lists have different lengths or if a mutation for any of these positions
                        already exists.
        """
        # Instead of checking *any* mutation, retrieve all mutation relationships between these proteins.
        # Here molecule1.mutation.relationship(molecule2) returns a list of mutation relationship instances.
        existing_mutations = molecule1.mutation.relationship(molecule2)

        if existing_mutations:
            raise ValueError(
                "A mutation relationship affecting one or more of these positions already exists between these proteins."
            )

        if (
            len(from_positions) != len(to_positions)
            or len(from_positions) != len(from_monomers)
            or len(from_positions) != len(to_monomers)
        ):
            raise ValueError("All input lists must have the same length.")

        for from_position, from_monomer in zip(from_positions, from_monomers):
            if molecule1.sequence[from_position] != from_monomer:
                raise ValueError(
                    f"Monomer '{from_monomer}' does not match the sequence {molecule1.accession_id} at position {from_position}"
                )

        for to_position, to_monomer in zip(to_positions, to_monomers):
            if molecule2.sequence[to_position] != to_monomer:
                raise ValueError(
                    f"Monomer '{to_monomer}' does not match the sequence {molecule2.accession_id} at position {to_position}"
                )

        molecule1.mutation.connect(
            molecule2,
            {
                "from_positions": from_positions,
                "to_positions": to_positions,
                "from_monomers": from_monomers,
                "to_monomers": to_monomers,
            },
        )

        return cls(
            from_positions=from_positions,
            to_positions=to_positions,
            from_monomers=from_monomers,
            to_monomers=to_monomers,
        )

    @property
    def label(self) -> str:
        """The label of the mutation."""
        return ",".join(
            f"{from_monomer}{from_position}{to_monomer}"
            for from_position, from_monomer, to_monomer in zip(
                list(self.from_positions),
                list(self.from_monomers),
                list(self.to_monomers),
            )
        )


class StandardNumberingRel(StructuredRel):  # type: ignore
    positions = ArrayProperty(StringProperty(), required=True)

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: Any,
        molecule2: Any,
        positions: list[str],
    ) -> "StandardNumberingRel":
        """Validates the positions and connects the two molecules."""
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
    def label(self) -> str:
        """The label of the standard numbering."""
        return f"{self.positions}"


class SiteRel(StructuredRel):  # type: ignore
    positions = ArrayProperty(IntegerProperty(), required=True)

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: Any,
        molecule2: Any,
        positions: list[int],
    ) -> "SiteRel":
        """Validates the positions and connects the two molecules."""
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
    def label(self) -> str:
        return f"{self.positions}"


class Site(StrictStructuredNode):
    site_id = UniqueIdProperty()
    name = StringProperty()
    annotation = StringProperty(
        choices=[(e.value, e.name) for e in Annotation], required=True
    )


class Region(StrictStructuredNode):
    region_id = UniqueIdProperty()
    name = StringProperty()
    annotation = StringProperty(
        choices=[(e.value, e.name) for e in Annotation], required=True
    )
    sequence_id = StringProperty()

    # Relationships
    has_mutation_region = RelationshipTo("Region", "MUTATION", model=Mutation)
    has_standard_numbering = RelationshipTo(
        "StandardNumbering", "HAS_STANDARD_NUMBERING", model=StandardNumberingRel
    )


class DNAProteinRel(StructuredRel):  # type: ignore
    """A relationship between a DNA and a protein."""

    start = IntegerProperty(required=True)
    end = IntegerProperty(required=True)

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: Any,
        molecule2: Any,
        start: int,
        end: int,
    ) -> "DNAProteinRel":
        """Validates the start and end positions and connects the two molecules."""
        molecule1.protein.connect(
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


class RegionRel(StructuredRel):  # type: ignore
    start = IntegerProperty(required=True)
    end = IntegerProperty(required=True)

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: Any,
        molecule2: Any,
        start: int,
        end: int,
    ) -> "RegionRel":
        """Validates the start and end positions and connects the two molecules."""
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
    def label(self) -> str:
        """The label of the region."""
        return f"{self.start}-{self.end}"


class Reaction(StrictStructuredNode):
    """
    A node representing a reaction.
    """

    rhea_id = StringProperty(unique_index=True, required=True)
    chebi_id = ArrayProperty(StringProperty())

    # Relationships
    substrate = RelationshipTo("Molecule", "SUBSTRATE")
    product = RelationshipTo("Molecule", "PRODUCT")

    @property
    def label(self) -> str:
        """The label of the reaction."""
        return f"{self.rhea_id}"


class Molecule(StrictStructuredNode):
    """
    A node representing a molecule in the database.
    """

    chebi_id = StringProperty(unique_index=True, required=True)
    rhea_compound_id = StringProperty()
    smiles = StringProperty()

    @classmethod
    def get_or_save(cls, **kwargs: Any) -> "Molecule":
        chebi_id = kwargs.get("chebi_id")
        smiles = kwargs.get("smiles")
        try:
            molecule = cast(Molecule, cls.nodes.get(chebi_id=chebi_id))
            return molecule
        except cls.DoesNotExist:
            try:
                molecule = cls(chebi_id=chebi_id, smiles=smiles)
                molecule.save()
                return molecule
            except Exception as e:
                print(f"Error during saving of the molecule: {e}")
                raise

    @property
    def label(self) -> str:
        """The label of the molecule."""
        return f"{self.chebi_id}"


class StandardNumbering(StrictStructuredNode):
    name = StringProperty(required=True, unique_index=True)
    definition = StringProperty(required=True)

    # Relationships
    sequences_protein = RelationshipTo(
        "Protein", "HAS_STANDARD_NUMBERING", model=StandardNumberingRel
    )


class PairwiseAlignmentResult(StructuredRel):  # type: ignore
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
        molecule1: Any,
        molecule2: Any,
        similarity: float,
        gaps: int,
        mismatches: int,
        score: int,
        query_aligned: str,
        target_aligned: str,
    ) -> "PairwiseAlignmentResult":
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
    def label(self) -> str:
        """The label of the GO annotation."""
        return str(self.term)


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
    rhea_id = ArrayProperty(StringProperty())
    chebi_id = ArrayProperty(StringProperty())
    embedding = ArrayProperty(
        FloatProperty(),
        vector_index=VectorIndex(dimensions=1280),
        index_type="hnsw",
        distance_metric="COSINE",
    )
    TBT = StringProperty()
    PCL = StringProperty()
    BHET = StringProperty()
    PET_powder = StringProperty()

    # Relationships
    organism = RelationshipTo("Organism", "ORIGINATES_FROM")
    site = RelationshipTo("Site", "HAS_SITE", model=SiteRel)
    region = RelationshipTo("Region", "HAS_REGION", model=RegionRel)
    go_annotation = RelationshipTo("GOAnnotation", "ASSOCIATED_WITH")
    reaction = RelationshipTo("Reaction", "HAS_REACTION")
    substrate = RelationshipTo("Molecule", "SUBSTRATE")
    product = RelationshipTo("Molecule", "PRODUCT")
    ontology_object = RelationshipTo("OntologyObject", "ASSOCIATED_WITH")
    mutation = RelationshipTo("Protein", "MUTATION", model=Mutation)
    pairwise_aligned = RelationshipTo(
        "Protein", "PAIRWISE_ALIGNED", model=PairwiseAlignmentResult
    )
    has_standard_numbering = RelationshipTo(
        "StandardNumbering", "HAS_STANDARD_NUMBERING", model=StandardNumberingRel
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
    protein = RelationshipTo("Protein", "ENCODES", model=DNAProteinRel)
    pairwise_aligned = RelationshipTo(
        "DNA", "PAIRWISE_ALIGNED", model=PairwiseAlignmentResult
    )
    has_standard_numbering = RelationshipTo(
        "StandardNumbering", "HAS_STANDARD_NUMBERING", model=StandardNumberingRel
    )


class CustomRealationship(StructuredRel):  # type: ignore
    """A custom relationship between two ontology objects."""

    name = StringProperty(required=True)
    description = StringProperty()

    @classmethod
    def validate_and_connect(
        cls,
        molecule1: Any,
        molecule2: Any,
        name: str,
        description: str,
    ) -> "CustomRealationship":
        molecule1.custom_relationships.connect(
            molecule2,
            {
                "name": name,
                "description": description,
            },
        )

        return cls(
            name=name,
            description=description,
        )

    @property
    def label(self) -> str:
        """The label of the custom relationship."""
        return str(self.name)


class OntologyObject(StrictStructuredNode):
    """A node representing an ontology object in the database."""

    name = StringProperty(required=True, unique_index=True)
    description = StringProperty()
    label = StringProperty()
    synonyms = ArrayProperty(StringProperty())

    # Relationships
    subclasses = RelationshipTo("OntologyObject", "SUBCLASS_OF")
    custom_relationships = RelationshipTo(
        "OntologyObject", "CUSTOM_RELATIONSHIP", model=CustomRealationship
    )
