from enum import Enum
from typing import List, Optional

from neomodel import (
    ArrayProperty,
    FloatProperty,
    IntegerProperty,
    RelationshipTo,
    StringProperty,
    StructuredNode,
    UniqueIdProperty,
    VectorIndex,
)


class Annotation(Enum):
    ACTIVE_SITE = "http://semanticscience.org/resource/SIO_010041"
    ALLOSTERIC_SITE = "http://semanticscience.org/resource/SIO_010050"
    ALPHAHELIX = "http://semanticscience.org/resource/SIO_010468"
    BETASTRAND = "http://semanticscience.org/resource/SIO_010469"
    BINDING_SITE = "http://semanticscience.org/resource/SIO_010040"
    CODING_SEQ = "http://semanticscience.org/resource/SIO_001276"
    DNA = "http://semanticscience.org/resource/SIO_010018"
    DOMAIN = "http://semanticscience.org/resource/SIO_001379"
    FAMILY = "http://semanticscience.org/resource/SIO_001380"
    MOTIVE = "http://semanticscience.org/resource/SIO_000131"
    PROTEIN = "http://semanticscience.org/resource/SIO_010015"


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


# Define the Site class
class Site(StructuredNode):
    uid = UniqueIdProperty()
    positions = ArrayProperty(IntegerProperty(), required=True)
    annotation = StringProperty(
        choices=[(e.value, e.name) for e in Annotation], required=True
    )


# Define the Region class
class Region(StructuredNode):
    region_id = UniqueIdProperty()
    start = IntegerProperty()
    end = IntegerProperty()
    annotation = StringProperty()


class ProteinRecord(StructuredNode):
    accession_id = StringProperty(unique_index=True, required=True)
    sequence = StringProperty(required=True)
    name = StringProperty()
    seq_length = IntegerProperty()
    mol_weight = FloatProperty()
    ec_number = StringProperty()
    nucleotide_id = StringProperty()
    locus_tag = StringProperty()
    structure_ids = ArrayProperty(StringProperty())
    go_terms = ArrayProperty(StringProperty())
    embedding = ArrayProperty(
        FloatProperty(),
        vector_index=VectorIndex(dimensions=1048, similarity_function="euclidean"),
    )

    # Relationships
    organism = RelationshipTo("Organism", "ORIGINATES_FROM")
    sites = RelationshipTo("Site", "HAS_SITE")
    regions = RelationshipTo("Region", "HAS_REGION")
    annotations = RelationshipTo("Annotation", "HAS_ANNOTATION")

    def add_to_sites(
        self,
        annotation: Annotation,
        name: Optional[str] = None,
        positions: List[int] = [],
        **kwargs,
    ):
        """Adds a new Site to the sites list"""
        params = {"annotation": annotation, "name": name, "positions": positions}
        if "id" in kwargs:
            params["id"] = kwargs["id"]
        site = Site(**params).save()
        self.sites.connect(site)
        return site

    def add_to_regions(
        self,
        region_id: str,
        annotation: Annotation,
        start: int,
        end: int,
        name: Optional[str] = None,
        **kwargs,
    ):
        """Adds a new Region to the regions list"""
        params = {
            "region_id": region_id,
            "annotation": annotation,
            "start": start,
            "end": end,
            "name": name,
        }
        region = Region(**params).save()
        self.regions.connect(region)
        return region

    def define_organism(self, taxonomy_id: int, **kwargs):
        """Defines the organism for the ProteinRecord"""
        params = {"taxonomy_id": taxonomy_id}
        params = params | kwargs

        organism = Organism(**params).save()
        self.organism.connect(organism)
        return organism


if __name__ == "__main__":
    # Create an instance of the ProteinRecord class
    from neomodel import config, db

    config.DATABASE_URL = "bolt://neo4j:none@localhost:7687"

    results, meta = db.cypher_query("RETURN 'Hello World' as message")
    print(results, meta)

    organism = Organism(taxonomy_id=960, name="ecol").save()

    prot = ProteinRecord(
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
