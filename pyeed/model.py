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


class Organism(StructuredNode):
    taxonomy_id = IntegerProperty(required=True, unique_index=True)
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


class Site(StructuredNode):
    site_id = UniqueIdProperty()
    positions = ArrayProperty(IntegerProperty(), required=True)
    annotation = StringProperty(
        choices=[(e.value, e.name) for e in Annotation], required=True
    )


class Region(StructuredNode):
    region_id = UniqueIdProperty()
    start = IntegerProperty(required=True)
    end = IntegerProperty(required=True)
    annotation = StringProperty(
        choices=[(e.value, e.name) for e in Annotation], required=True
    )


class Protein(StructuredNode):
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


class DNA(StructuredNode):
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
