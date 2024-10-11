from enum import Enum

from neomodel import (
    ArrayProperty,
    FloatProperty,
    IntegerProperty,
    RelationshipFrom,
    RelationshipTo,
    StringProperty,
    StructuredNode,
    UniqueIdProperty,
    VectorIndex,
)


class StrictStructuredNode(StructuredNode):
    """A StructuredNode subclass that raises an error if an invalid property is provided."""

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
        for field, prop in self.__dict__.items():
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
    domain = StringProperty()
    kingdom = StringProperty()
    phylum = StringProperty()
    tax_class = StringProperty()
    order = StringProperty()
    family = StringProperty()
    genus = StringProperty()
    species = StringProperty()

    # Relationships
    protein = RelationshipFrom("Protein", "ORIGINATES_FROM")
    dna = RelationshipFrom("DNA", "ORIGINATES_FROM")

    @classmethod
    def add_or_skip(cls, **kwargs):
        """Add an organism if it does not already exist."""
        taxonomy_id = kwargs.get("taxonomy_id")
        organism = cls.nodes.get_or_none(taxonomy_id=taxonomy_id)
        if organism is None:
            organism = cls(**kwargs).save()

        return organism


class Site(StrictStructuredNode):
    site_id = UniqueIdProperty()
    name = StringProperty()
    positions = ArrayProperty(IntegerProperty(), required=True)
    annotation = StringProperty(
        choices=[(e.value, e.name) for e in Annotation], required=True
    )


class Region(StrictStructuredNode):
    region_id = UniqueIdProperty()
    start = IntegerProperty(required=True)
    end = IntegerProperty(required=True)
    annotation = StringProperty(
        choices=[(e.value, e.name) for e in Annotation], required=True
    )


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
        vector_index=VectorIndex(dimensions=1048),
    )

    # Relationships
    organism = RelationshipTo("Organism", "ORIGINATES_FROM")
    sites = RelationshipTo("Site", "HAS_SITE")
    regions = RelationshipTo("Region", "HAS_REGION")


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
    sites = RelationshipTo("Site", "HAS_SITE")
    regions = RelationshipTo("Region", "HAS_REGION")


if __name__ == "__main__":
    # Create an instance of the Protein class
    from neomodel import config, db

    config.DATABASE_URL = "bolt://neo4j:none@localhost:7687"

    results, meta = db.cypher_query("RETURN 'Hello World' as message")
    print(results, meta)

    organism = Organism(taxonomy_id=960, name="ecol").save()

    prot = Protein(
        accession_id="P12345",
        sequence="MSEQ1234",
        name="Protein A",
        mol_weight=123.45,
    ).save()

    prot.organism.connect(organism)

    prot.add_to_regions(
        region_id="R123",
        annotation="Region",
        start=1,
        end=10,
    )
    prot.add_to_regions(
        region_id="R124",
        annotation="Region",
        start=11,
        end=20,
    )
