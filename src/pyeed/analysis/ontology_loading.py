from Bio.Align import Alignment as Alignment
from pyeed.dbconnect import DatabaseConnector
from pyeed.main import Pyeed
from rdflib import Graph, RDF, RDFS, OWL, Namespace



class OntologyAdapter():
    """
    Adapter class to load ontology files into the database.
    """

    def import_ontology_file_in_db(self, file_path: str, db: DatabaseConnector):
        """
        Imports an ontology file into the database.

        :param file_path: The path to the ontology file.
        :param db: The database connector

        :return: None
        """
        
        # Load the OWL file
        g = Graph()
        g.parse(file_path)

        # Create a namespace for the ontology
        IAO_NS = Namespace("http://purl.obolibrary.org/obo/IAO_")
        OBOINOWL_NS = Namespace("http://www.geneontology.org/formats/oboInOwl#")

        # create a dictonary of the labels
        dicts_labels = {}
        for s, p, o in g.triples((None, RDFS.label, None)):
            dicts_labels[str(s)] = str(o)
    
        # Iterate over the classes in the OWL file
        for s, p, o in g.triples((None, RDF.type, OWL.Class)):
            class_name = str(s)
            db.execute_write("CREATE (c:OntologyObject {name: $name})", parameters = {"name": class_name})

            # add discreption, example in CARD: <obo:IAO_0000115>eccC5 is a.....</obo:IAO_0000115>
            for _, _, desc in g.triples((s, IAO_NS['0000115'], None)):
                description = str(desc)
                db.execute_write("""
                    MATCH (c:OntologyObject {name: $name})
                    SET c.description = $description
                """, parameters = {"name": class_name, "description": description})

            # add the label to the class
            db.execute_write("""
                MATCH (c:OntologyObject {name: $name})
                SET c.label = $label
            """, parameters = {"name": class_name, "label": dicts_labels[class_name]})

            # add the synonyms to the class
            # <oboInOwl:hasExactSynonym>Mtub_eccC5_FLO</oboInOwl:hasExactSynonym>
            for _, _, syn in g.triples((s, OBOINOWL_NS.hasExactSynonym, None)):
                synonym = str(syn)
                db.execute_write("""
                    MATCH (c:OntologyObject {name: $name})
                    SET c.synonym = $synonym
                """, parameters = {"name": class_name, "synonym": synonym})


        # Create relationships (subclasses, properties)
        for s, p, o in g.triples((None, RDFS.subClassOf, None)):
            if (o, RDF.type, OWL.Class) in g:
                subclass = str(s)
                superclass = str(o)
                db.execute_write("""
                    MATCH (sub:OntologyObject {name: $subclass}), (super:OntologyObject {name: $superclass})
                    CREATE (sub)-[:SUBCLASS_OF]->(super)
                """, parameters = {"subclass": subclass, "superclass": superclass})

            # handels the case where the subclass is a restriction, RO_ (in CARD)
            elif (o, RDF.type, OWL.Restriction) in g:
                on_property = None
                some_values_from = None
                
                # Extract onProperty
                for _, _, prop in g.triples((o, OWL.onProperty, None)):
                    on_property = str(prop)
                
                # Extract someValuesFrom
                for _, _, value in g.triples((o, OWL.someValuesFrom, None)):
                    some_values_from = str(value)
                
                if on_property and some_values_from:
                    # create a realtionship of type CustomRealationship with the name on_property and the description which can be checked in the dict
                    # link is between the subclass and the some_values_from
                    db.execute_write("""
                        MATCH (sub:OntologyObject {name: $subclass}), (super:OntologyObject {name: $some_values_from})
                        CREATE (sub)-[:CustomRelationship {name: $on_property, description: $description}]->(super)
                    """, parameters = {"subclass": subclass, "some_values_from": some_values_from, "on_property": on_property, "description": dicts_labels[on_property]})





if __name__ == "__main__":

    uri = "bolt://localhost:7687"
    username = "neo4j"
    password = "12345678"

    file_path = "/home/nab/Niklas/TEM-lactamase/CARD_Data_Ontologies/aro.owl"

    eedb = Pyeed(uri, user=username, password=password)
    eedb.db.wipe_database()
    eedb.db.remove_db_constraints(user=username, password=password)

    eedb.db.initialize_db_constraints(user=username, password=password)


    db = eedb.db
    ontology_adapter = OntologyAdapter()

    ontology_adapter.import_ontology_file_in_db(file_path, db)
    