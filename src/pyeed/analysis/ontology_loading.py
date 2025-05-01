from typing import Dict

from pyeed.dbconnect import DatabaseConnector
from rdflib import OWL, RDF, RDFS, Graph, Namespace, URIRef


class OntologyAdapter:
    """
    Adapter class to load ontology files into the database.
    Handles parsing and importing of OWL ontology files into a Neo4j database.
    """

    def import_ontology_file_in_db(self, file_path: str, db: DatabaseConnector) -> None:
        """
        Imports an ontology file into the Neo4j database.

        This method parses an OWL file and creates corresponding nodes and relationships
        in the Neo4j database. It handles class definitions, labels, descriptions,
        synonyms, and relationships including subclass relationships and custom
        relationships defined through OWL restrictions.

        Args:
            file_path: Path to the OWL ontology file to import
            db: DatabaseConnector instance for Neo4j database operations

        Raises:
            FileNotFoundError: If the ontology file cannot be found
            rdflib.exceptions.ParserError: If the OWL file cannot be parsed
        """
        # Load the OWL file
        g = Graph()
        g.parse(file_path)

        # Create namespaces for the ontology
        IAO_NS = Namespace("http://purl.obolibrary.org/obo/IAO_")
        OBOINOWL_NS = Namespace("http://www.geneontology.org/formats/oboInOwl#")

        # Create a dictionary of the labels
        dicts_labels: Dict[str, str] = {}
        for s, _, o in g.triples((None, RDFS.label, None)):
            dicts_labels[str(s)] = str(o)

        # Process OWL classes
        self._process_owl_classes(g, db, dicts_labels, IAO_NS, OBOINOWL_NS)

        # Process relationships
        self._process_relationships(g, db, dicts_labels)

    def _process_owl_classes(
        self,
        g: Graph,
        db: DatabaseConnector,
        dicts_labels: Dict[str, str],
        IAO_NS: Namespace,
        OBOINOWL_NS: Namespace,
    ) -> None:
        """Process OWL classes and create corresponding database nodes."""
        for s_node, p_node, o_node in g.triples((None, RDF.type, OWL.Class)):
            class_name = str(s_node)
            # Create node for the class
            db.execute_write(
                "CREATE (c:OntologyObject {name: $name})",
                parameters={"name": class_name},
            )

            # Add description if available
            for _, _, desc in g.triples((s_node, IAO_NS["0000115"], None)):
                description = str(desc)
                db.execute_write(
                    """
                    MATCH (c:OntologyObject {name: $name})
                    SET c.description = $description
                    """,
                    parameters={"name": class_name, "description": description},
                )

            # Add label if available
            if class_name in dicts_labels:
                db.execute_write(
                    """
                    MATCH (c:OntologyObject {name: $name})
                    SET c.label = $label
                    """,
                    parameters={"name": class_name, "label": dicts_labels[class_name]},
                )

            # Add synonyms
            for _, _, syn in g.triples((s_node, OBOINOWL_NS.hasExactSynonym, None)):
                synonym = str(syn)
                db.execute_write(
                    """
                    MATCH (c:OntologyObject {name: $name})
                    SET c.synonym = $synonym
                    """,
                    parameters={"name": class_name, "synonym": synonym},
                )

    def _process_relationships(
        self,
        g: Graph,
        db: DatabaseConnector,
        dicts_labels: Dict[str, str],
    ) -> None:
        """Process OWL relationships and create corresponding database relationships."""
        for s, p, o in g.triples((None, RDFS.subClassOf, None)):
            subclass = str(s)

            if (o, RDF.type, OWL.Class) in g:
                # Handle direct subclass relationships
                superclass = str(o)
                db.execute_write(
                    """
                    MATCH (sub:OntologyObject {name: $subclass}),
                          (super:OntologyObject {name: $superclass})
                    CREATE (sub)-[:SUBCLASS_OF]->(super)
                    """,
                    parameters={"subclass": subclass, "superclass": superclass},
                )

            elif (o, RDF.type, OWL.Restriction) in g:
                # Handle OWL restrictions (e.g., RO_ in CARD)
                self._process_restriction(g, str(o), subclass, db, dicts_labels)

    def _process_restriction(
        self,
        g: Graph,
        restriction_node: str,
        subclass: str,
        db: DatabaseConnector,
        dicts_labels: Dict[str, str],
    ) -> None:
        """Process OWL restrictions and create custom relationships."""
        on_property = None
        some_values_from = None

        # Convert restriction_node string to RDFLib URIRef
        restriction = URIRef(restriction_node)

        # Extract onProperty
        for _, _, prop in g.triples((restriction, OWL.onProperty, None)):
            on_property = str(prop)

        # Extract someValuesFrom
        for _, _, value in g.triples((restriction, OWL.someValuesFrom, None)):
            some_values_from = str(value)

        if on_property and some_values_from:
            query_params = {
                "subclass": subclass,
                "some_values_from": some_values_from,
                "on_property": on_property,
                "description": dicts_labels.get(on_property, ""),
            }

            query = """
                MATCH (sub:OntologyObject {name: $subclass}),
                      (super:OntologyObject {name: $some_values_from})
                CREATE (sub)-[:CustomRelationship {
                    name: $on_property,
                    description: $description
                }]->(super)
            """

            db.execute_write(query, parameters=query_params)
