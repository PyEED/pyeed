import asyncio
import logging
from enum import Enum
from typing import Any, Dict, List, Optional, Set, TypeVar

from neo4j import GraphDatabase
from pydantic import BaseModel

from .model import (
    BaseNode,
    Embedding,
    GOAnnotation,
    NodeHint,
    Organism,
    Protein,
    Reaction,
    SequenceAnnotation,
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# add file handler if not already present
if not logger.handlers:
    fh = logging.FileHandler("neo4j_integration.log", mode="a")  # "a" for append
    fh.setLevel(logging.DEBUG)

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    fh.setFormatter(formatter)

    logger.addHandler(fh)

T = TypeVar("T", bound=BaseNode)

_PRIMITIVES = (str, int, float, bool)


def extract_node_hints(model: type[BaseModel]) -> Dict[str, NodeHint]:
    """Extract NodeHint metadata from model fields."""
    hints: Dict[str, NodeHint] = {}
    logger.debug(f"Extracting NodeHints from model: {model.__name__}")

    for field_name, field_info in model.model_fields.items():
        for meta in getattr(field_info, "metadata", ()):
            if isinstance(meta, NodeHint):
                hints[field_name] = meta
                logger.debug(f"Found NodeHint for field '{field_name}': {meta}")
                break

    logger.debug(f"Extracted {len(hints)} NodeHints: {hints}")
    return hints


# EdgeHint functionality removed - using manual edge creation instead


def get_model_class_by_name(class_name: str) -> Optional[type[BaseNode]]:
    """Get a model class by name without hardcoded imports."""
    logger.debug(f"Looking for model class: {class_name}")
    try:
        # Dynamic import to get the class
        import sys

        base_module = sys.modules.get("pyeed.model")
        if base_module:
            model_class = getattr(base_module, class_name, None)
            if model_class:
                logger.debug(f"Found model class: {class_name}")
            else:
                logger.debug(f"Model class not found: {class_name}")
            return model_class
        else:
            logger.debug("pyeed.model module not found")
    except Exception as e:
        logger.debug(f"Error getting model class {class_name}: {e}")
    return None


class Neo4jNode:
    """Represents a Neo4j node with its properties and labels."""

    def __init__(
        self,
        labels: Set[str],
        properties: Dict[str, Any],
        node_hints: Dict[str, Any],
        unique_id: str,
    ):
        self.labels = labels
        self.properties = properties
        self.node_hints = node_hints
        self.unique_id: str = unique_id

    @property
    def id(self) -> str:
        """Return the unique identifier for this node."""
        return self.unique_id

    def __repr__(self) -> str:
        return (
            f"Neo4jNode(\n"
            f"  labels={self.labels},\n"
            f"  properties={self.properties},\n"
            f"  node_hints={self.node_hints},\n"
            f"  unique_id={self.unique_id}\n"
            f")"
        )


class Neo4jRelationship:
    """Represents a Neo4j relationship between two nodes."""

    def __init__(
        self,
        start_node_id: str,
        end_node_id: str,
        rel_type: str,
        properties: Dict[str, Any],
    ):
        self.start_node_id = start_node_id
        self.end_node_id = end_node_id
        self.rel_type = rel_type
        self.properties = properties

    def __repr__(self) -> str:
        return f"Neo4jRelationship({self.start_node_id})-[:{self.rel_type}]->({self.end_node_id})"


class PyeedNeo4jMapper:
    """Generic mapper for Pyeed models using NodeHint and EdgeHint metadata."""

    def __init__(self, driver: Any):
        logger.info("Initializing PyeedNeo4jMapper")
        self.driver = driver
        logger.info("PyeedNeo4jMapper initialized successfully")

    def _create_node_properties(self, model: BaseNode) -> Dict[str, Any]:
        """Create node properties from a Pydantic model, flattening custom fields."""
        properties: Dict[str, Any] = {}
        logger.debug(f"Creating node properties for model: {model.__class__.__name__}")

        for field_name, field_value in model.model_dump().items():
            if field_value is None:
                logger.debug(f"Skipping field '{field_name}': is None")
                continue

            # Skip private fields except for _id which we need for relationships
            if field_name.startswith("_") and field_name != "_id":
                logger.debug(f"Skipping private field '{field_name}'")
                continue

            # Skip dicts entirely
            if isinstance(field_value, dict):
                logger.debug(f"Skipping dict field '{field_name}': {type(field_value)}")
                continue

            # Only include non-empty lists of primitives
            if isinstance(field_value, list):
                if not field_value:
                    logger.debug(f"Skipping empty list field '{field_name}'")
                    continue
                if all(isinstance(el, _PRIMITIVES) or el is None for el in field_value):
                    properties[field_name] = field_value
                    logger.debug(f"Added list field '{field_name}': {field_value}")
                else:
                    logger.debug(
                        f"Skipping non-primitive list field '{field_name}': {field_value}"
                    )
                continue

            # If enum, use value
            if isinstance(field_value, Enum):
                properties[field_name] = field_value.value
                logger.debug(f"Added enum field '{field_name}': {field_value.value}")
                continue

            # Only include primitives
            if isinstance(field_value, _PRIMITIVES):
                properties[field_name] = field_value
                logger.debug(f"Added primitive field '{field_name}': {field_value}")
                continue

            logger.debug(
                f"Skipping field '{field_name}': {type(field_value)} - {field_value}"
            )

        # Merge in custom fields
        if model.custom:
            logger.debug(f"Adding custom fields: {model.custom}")
            for key, value in model.custom.items():
                if key in properties:
                    msg = f"Key {key} already exists in properties"
                    logger.error(msg)
                    raise ValueError(msg)
                properties[key] = value
                logger.debug(f"Added custom field '{key}': {value}")

        # Ensure _id is always included for relationship creation
        if hasattr(model, "_id") and "_id" not in properties:
            properties["_id"] = model._id
            logger.debug(f"Added _id field: {model._id}")

        logger.debug(f"Final properties for {model.__class__.__name__}: {properties}")
        return properties

    def create_node(self, model: BaseNode) -> Neo4jNode:
        """Create a Neo4j node from any BaseNode model."""
        logger.debug(f"Creating Neo4j node for model: {model.__class__.__name__}")

        # Get the unique ID directly from the model
        unique_id = model.get_unique_id()
        logger.debug(f"Model {model.__class__.__name__} unique_id: {unique_id}")
        logger.debug(f"Model {model.__class__.__name__} _id: {model._id}")

        node_hints = extract_node_hints(type(model))
        logger.debug(f"Node hints for {model.__class__.__name__}: {node_hints}")

        # Create main node
        labels = {model.__class__.__name__}
        properties = self._create_node_properties(model)
        node = Neo4jNode(labels, properties, node_hints, unique_id)

        logger.debug(f"Created Neo4jNode: {node}")
        return node

    def create_relationships_manual(
        self, source_model: BaseNode, target_models: List[BaseNode], rel_type: str
    ) -> List[Neo4jRelationship]:
        """Create relationships manually without EdgeHint metadata."""
        logger.debug(
            f"Creating manual relationships from {source_model.__class__.__name__} to {len(target_models)} targets"
        )
        logger.debug(f"Relationship type: {rel_type}")

        relationships = []

        # Get source node unique ID directly
        source_id = source_model.get_unique_id()
        logger.debug(
            f"Source model {source_model.__class__.__name__} unique_id: {source_id}"
        )

        for i, target_model in enumerate(target_models):
            # Get target node unique ID directly
            target_id = target_model.get_unique_id()
            logger.debug(
                f"Target model {i+1} {target_model.__class__.__name__} unique_id: {target_id}"
            )

            # Create simple relationship properties
            rel_properties = {
                "source_id": source_id,
                "target_id": target_id,
            }
            logger.debug(f"Relationship properties: {rel_properties}")

            rel = Neo4jRelationship(
                start_node_id=source_id,
                end_node_id=target_id,
                rel_type=rel_type.upper(),
                properties=rel_properties,
            )
            logger.debug(f"Created relationship: {rel}")
            relationships.append(rel)

        logger.debug(f"Created {len(relationships)} relationships")
        return relationships

    def map_model_to_neo4j(self, model: BaseNode) -> Dict[str, Any]:
        """Map any BaseNode model to Neo4j nodes and relationships using manual mapping."""
        logger.debug(f"Mapping model to Neo4j: {model.__class__.__name__}")
        result: Dict[str, Any] = {"nodes": [], "relationships": []}

        # Create main model node
        logger.debug("Creating main model node...")
        main_node = self.create_node(model)
        result["nodes"].append(main_node)
        logger.debug(f"Added main node: {main_node}")

        # Manual relationship mapping based on known model structure
        if hasattr(model, "annotations") and model.annotations:
            logger.debug(f"Processing {len(model.annotations)} sequence annotations")
            for annotation in model.annotations:
                annotation_node = self.create_node(annotation)
                result["nodes"].append(annotation_node)
                logger.debug(f"Added annotation node: {annotation_node}")

            # Create relationships from protein to annotations
            relationships = self.create_relationships_manual(
                model, model.annotations, "HAS_ANNOTATION"
            )
            result["relationships"].extend(relationships)
            logger.debug(f"Added {len(relationships)} annotation relationships")

        if hasattr(model, "go_terms") and model.go_terms:
            logger.debug(f"Processing {len(model.go_terms)} GO terms")
            for go_term in model.go_terms:
                go_node = self.create_node(go_term)
                result["nodes"].append(go_node)
                logger.debug(f"Added GO term node: {go_node}")

            # Create relationships from protein to GO terms
            relationships = self.create_relationships_manual(
                model, model.go_terms, "HAS_GO_ANNOTATION"
            )
            result["relationships"].extend(relationships)
            logger.debug(f"Added {len(relationships)} GO term relationships")

        if hasattr(model, "organisms") and model.organisms:
            logger.debug(f"Processing {len(model.organisms)} organisms")
            for organism in model.organisms:
                organism_node = self.create_node(organism)
                result["nodes"].append(organism_node)
                logger.debug(f"Added organism node: {organism_node}")

            # Create relationships from protein to organisms
            relationships = self.create_relationships_manual(
                model, model.organisms, "ORIGINATES_FROM"
            )
            result["relationships"].extend(relationships)
            logger.debug(f"Added {len(relationships)} organism relationships")

        if hasattr(model, "embeddings") and model.embeddings:
            logger.debug(f"Processing {len(model.embeddings)} embeddings")
            for embedding in model.embeddings.values():
                embedding_node = self.create_node(embedding)
                result["nodes"].append(embedding_node)
                logger.debug(f"Added embedding node: {embedding_node}")

            # Create relationships from protein to embeddings
            relationships = self.create_relationships_manual(
                model, list(model.embeddings.values()), "HAS_EMBEDDING"
            )
            result["relationships"].extend(relationships)
            logger.debug(f"Added {len(relationships)} embedding relationships")

        if hasattr(model, "reactions") and model.reactions:
            logger.debug(f"Processing {len(model.reactions)} reactions")
            for reaction in model.reactions:
                reaction_node = self.create_node(reaction)
                result["nodes"].append(reaction_node)
                logger.debug(f"Added reaction node: {reaction_node}")

            # Create relationships from protein to reactions
            relationships = self.create_relationships_manual(
                model, model.reactions, "HAS_REACTION"
            )
            result["relationships"].extend(relationships)
            logger.debug(f"Added {len(relationships)} reaction relationships")

        logger.debug(
            f"Final mapping result: {len(result['nodes'])} nodes, {len(result['relationships'])} relationships"
        )
        return result

    def create_constraints_and_indexes(self) -> None:
        """Create Neo4j constraints and indexes based on NodeHint metadata."""
        logger.info("Creating Neo4j constraints and indexes")
        with self.driver.session() as session:
            # Dynamically discover all BaseNode subclasses without hardcoding imports
            import inspect
            import sys

            # Get the module where BaseNode is defined
            base_module = sys.modules.get("pyeed.model")
            if not base_module:
                logger.warning(
                    "Could not find pyeed.model module for dynamic discovery"
                )
                return

            # Find all classes that inherit from BaseNode
            model_classes = []
            for name, obj in inspect.getmembers(base_module):
                if (
                    inspect.isclass(obj)
                    and issubclass(obj, BaseNode)
                    and obj != BaseNode
                ):
                    model_classes.append(obj)
                    logger.debug(f"Found BaseNode subclass: {name}")

            logger.info(
                f"Discovered {len(model_classes)} BaseNode subclasses: {[cls.__name__ for cls in model_classes]}"
            )

            for model_class in model_classes:
                logger.debug(f"Processing model class: {model_class.__name__}")
                # Check each field for unique constraints
                node_hints = extract_node_hints(model_class)
                logger.debug(f"Node hints for {model_class.__name__}: {node_hints}")

                for field_name, hint in node_hints.items():
                    logger.debug(f"Processing hint for field {field_name}: {hint}")

                    if hint.unique:
                        label = model_class.__name__
                        logger.debug(
                            f"Creating unique constraint for {label}.{field_name}"
                        )
                        try:
                            cypher = f"CREATE CONSTRAINT {label.lower()}_{field_name}_unique IF NOT EXISTS FOR (n:{label}) REQUIRE n.{field_name} IS UNIQUE"
                            logger.debug(f"Constraint cypher: {cypher}")
                            session.run(cypher)
                            logger.info(
                                f"Created unique constraint on {label}.{field_name}"
                            )
                        except Exception as e:
                            logger.warning(
                                f"Could not create constraint on {label}.{field_name}: {e}"
                            )

                    if hint.vector_index:
                        label = model_class.__name__
                        logger.debug(f"Creating vector index for {label}.{field_name}")
                        try:
                            cypher = f"CREATE VECTOR INDEX {label.lower()}_{field_name}_vector IF NOT EXISTS FOR (n:{label}) ON (n.{field_name}) OPTIONS {{indexConfig: {{`vector.dimensions`: 1280, `vector.similarity_function`: 'cosine'}}}}"
                            logger.debug(f"Vector index cypher: {cypher}")
                            session.run(cypher)
                            logger.info(f"Created vector index on {label}.{field_name}")
                        except Exception as e:
                            logger.warning(
                                f"Could not create vector index on {label}.{field_name}: {e}"
                            )


class PyeedNeo4jWriter:
    """Generic writer for Pyeed models to Neo4j database."""

    def __init__(self, uri: str, user: str, password: str):
        logger.info(f"Initializing PyeedNeo4jWriter with URI: {uri}")
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
        self.mapper = PyeedNeo4jMapper(self.driver)
        logger.info("PyeedNeo4jWriter initialized successfully")

    def close(self) -> None:
        """Close the Neo4j driver connection."""
        logger.info("Closing Neo4j driver connection")
        self.driver.close()
        logger.info("Neo4j driver connection closed")

    def _create_node(self, model: BaseNode) -> Dict[str, Any]:
        """Helper method to create a single node in Neo4j."""
        logger.debug(f"Creating node for {model.__class__.__name__}")

        with self.driver.session() as session:
            neo4j_node = self.mapper.create_node(model)

            # Always try to match first, then create if needed using MERGE
            labels_str = ":".join(neo4j_node.labels)
            unique_field = model.get_unique_field_name()
            unique_value = model.get_unique_id()

            # Use MERGE to match existing or create new node
            cypher = f"""
            MERGE (n:{labels_str} {{{unique_field}: $unique_value}})
            ON CREATE SET n = $properties
            ON MATCH SET n += $properties
            RETURN n
            """

            logger.debug(f"Cypher: {cypher}")
            logger.debug(f"Properties: {neo4j_node.properties}")
            logger.debug(f"Unique field: {unique_field}, value: {unique_value}")

            result = session.run(
                cypher,
                {"unique_value": unique_value, "properties": neo4j_node.properties},
            )
            created_node = result.single()

            if not created_node:
                logger.error(f"Failed to merge node: {neo4j_node}")
                return {"success": False, "error": "Failed to merge node"}

            neo4j_id = str(created_node["n"].element_id)
            logger.info(f"Merged node {neo4j_node.labels} with Neo4j ID {neo4j_id}")

            return {"success": True, "neo4j_id": neo4j_id}

    def add_protein(self, protein: Protein) -> Dict[str, Any]:
        """Add a protein with all its nested objects to Neo4j."""

        if not isinstance(protein, Protein):
            return {"success": False, "error": "Object must be a Protein instance"}

        logger.info(f"Adding protein {protein.accession_id} to Neo4j")

        try:
            with self.driver.session() as session:
                # Create the main protein node
                protein_result = self._create_node(protein)
                if not protein_result["success"]:
                    return protein_result

                total_nodes = 1
                total_relationships = 0

                # Add all organisms
                for organism in protein.organisms:
                    org_result = self.add_organism(protein.accession_id, organism)
                    if org_result["success"]:
                        total_nodes += 1
                        total_relationships += 1

                # Add all sequence annotations
                for annotation in protein.annotations:
                    ann_result = self.add_sequence_annotation(
                        protein.accession_id, annotation
                    )
                    if ann_result["success"]:
                        total_nodes += 1
                        total_relationships += 1

                # Add all GO annotations
                for go_term in protein.go_terms:
                    go_result = self.add_go_annotation(protein.accession_id, go_term)
                    if go_result["success"]:
                        total_nodes += 1
                        total_relationships += 1

                # Add all embeddings
                for embedding in protein.embeddings.values():
                    emb_result = self.add_embedding(protein.accession_id, embedding)
                    if emb_result["success"]:
                        total_nodes += 1
                        total_relationships += 1

                # Add all reactions
                for reaction in protein.reactions:
                    react_result = self.add_reaction(protein.accession_id, reaction)
                    if react_result["success"]:
                        total_nodes += 1
                        total_relationships += 1

                logger.info(
                    f"Successfully added protein {protein.accession_id} with {total_nodes} nodes and {total_relationships} relationships"
                )
                return {
                    "success": True,
                    "nodes_created": total_nodes,
                    "relationships_created": total_relationships,
                }

        except Exception as e:
            logger.error(f"Error adding protein {protein.accession_id} to Neo4j: {e}")
            logger.exception("Full traceback:")
            return {"success": False, "error": str(e)}

    def add_organism(
        self, protein_accession_id: str, organism: Organism
    ) -> Dict[str, Any]:
        """Add an organism to Neo4j and link it to a protein."""
        from .model import Organism

        if not isinstance(organism, Organism):
            return {"success": False, "error": "Object must be an Organism instance"}

        logger.info(
            f"Adding organism {organism.tax_id} and linking to protein {protein_accession_id}"
        )

        try:
            with self.driver.session() as session:
                # Create the organism node
                org_result = self._create_node(organism)
                if not org_result["success"]:
                    return org_result

                # Create relationship to protein
                rel_cypher = """
                MATCH (p:Protein {accession_id: $protein_id})
                MATCH (o:Organism {tax_id: $organism_id})
                MERGE (p)-[r:ORIGINATES_FROM]->(o)
                RETURN r
                """

                session.run(
                    rel_cypher,
                    {
                        "protein_id": protein_accession_id,
                        "organism_id": organism.tax_id,
                    },
                )

                logger.info(
                    f"Created relationship between protein {protein_accession_id} and organism {organism.tax_id}"
                )
                return {"success": True, "relationship_created": True}

        except Exception as e:
            logger.error(f"Error adding organism: {e}")
            return {"success": False, "error": str(e)}

    def add_sequence_annotation(
        self, protein_accession_id: str, annotation: SequenceAnnotation
    ) -> Dict[str, Any]:
        """Add a sequence annotation to Neo4j and link it to a protein."""

        if not isinstance(annotation, SequenceAnnotation):
            return {
                "success": False,
                "error": "Object must be a SequenceAnnotation instance",
            }

        logger.info(
            f"Adding sequence annotation and linking to protein {protein_accession_id}"
        )

        try:
            with self.driver.session() as session:
                # Create the annotation node
                ann_result = self._create_node(annotation)
                if not ann_result["success"]:
                    return ann_result

                # Create relationship to protein
                rel_cypher = """
                MATCH (p:Protein {accession_id: $protein_id})
                MATCH (a:SequenceAnnotation) WHERE elementId(a) = $annotation_id
                MERGE (p)-[r:HAS_ANNOTATION]->(a)
                RETURN r
                """

                session.run(
                    rel_cypher,
                    {
                        "protein_id": protein_accession_id,
                        "annotation_id": ann_result["neo4j_id"],
                    },
                )

                logger.info(
                    f"Created relationship between protein {protein_accession_id} and sequence annotation"
                )
                return {"success": True, "relationship_created": True}

        except Exception as e:
            logger.error(f"Error adding sequence annotation: {e}")
            return {"success": False, "error": str(e)}

    def add_go_annotation(
        self, protein_accession_id: str, go_annotation: GOAnnotation
    ) -> Dict[str, Any]:
        """Add a GO annotation to Neo4j and link it to a protein."""

        if not isinstance(go_annotation, GOAnnotation):
            return {"success": False, "error": "Object must be a GOAnnotation instance"}

        logger.info(
            f"Adding GO annotation {go_annotation.go_id} and linking to protein {protein_accession_id}"
        )

        try:
            with self.driver.session() as session:
                # Create the GO annotation node
                go_result = self._create_node(go_annotation)
                if not go_result["success"]:
                    return go_result

                # Create relationship to protein
                rel_cypher = """
                MATCH (p:Protein {accession_id: $protein_id})
                MATCH (g:GOAnnotation {go_id: $go_id})
                MERGE (p)-[r:HAS_GO_ANNOTATION]->(g)
                RETURN r
                """

                session.run(
                    rel_cypher,
                    {"protein_id": protein_accession_id, "go_id": go_annotation.go_id},
                )

                logger.info(
                    f"Created relationship between protein {protein_accession_id} and GO annotation {go_annotation.go_id}"
                )
                return {"success": True, "relationship_created": True}

        except Exception as e:
            logger.error(f"Error adding GO annotation: {e}")
            return {"success": False, "error": str(e)}

    def add_embedding(
        self, protein_accession_id: str, embedding: Embedding
    ) -> Dict[str, Any]:
        """Add an embedding to Neo4j and link it to a protein."""

        if not isinstance(embedding, Embedding):
            return {"success": False, "error": "Object must be an Embedding instance"}

        logger.info(
            f"Adding embedding {embedding.description} and linking to protein {protein_accession_id}"
        )

        try:
            with self.driver.session() as session:
                # Create the embedding node
                emb_result = self._create_node(embedding)
                if not emb_result["success"]:
                    return emb_result

                # Create relationship to protein
                rel_cypher = """
                MATCH (p:Protein {accession_id: $protein_id})
                MATCH (e:Embedding) WHERE elementId(e) = $embedding_id
                MERGE (p)-[r:HAS_EMBEDDING]->(e)
                RETURN r
                """

                session.run(
                    rel_cypher,
                    {
                        "protein_id": protein_accession_id,
                        "embedding_id": emb_result["neo4j_id"],
                    },
                )

                logger.info(
                    f"Created relationship between protein {protein_accession_id} and embedding {embedding.description}"
                )
                return {"success": True, "relationship_created": True}

        except Exception as e:
            logger.error(f"Error adding embedding: {e}")
            return {"success": False, "error": str(e)}

    def add_reaction(
        self, protein_accession_id: str, reaction: Reaction
    ) -> Dict[str, Any]:
        """Add a reaction to Neo4j and link it to a protein."""

        if not isinstance(reaction, Reaction):
            return {"success": False, "error": "Object must be a Reaction instance"}

        logger.info(
            f"Adding reaction {reaction.rhea_id} and linking to protein {protein_accession_id}"
        )

        try:
            with self.driver.session() as session:
                # Create the reaction node
                react_result = self._create_node(reaction)
                if not react_result["success"]:
                    return react_result

                # Create relationship to protein
                rel_cypher = """
                MATCH (p:Protein {accession_id: $protein_id})
                MATCH (r:Reaction {rhea_id: $rhea_id})
                MERGE (p)-[rel:HAS_REACTION]->(r)
                RETURN rel
                """

                session.run(
                    rel_cypher,
                    {"protein_id": protein_accession_id, "rhea_id": reaction.rhea_id},
                )

                logger.info(
                    f"Created relationship between protein {protein_accession_id} and reaction {reaction.rhea_id}"
                )
                return {"success": True, "relationship_created": True}

        except Exception as e:
            logger.error(f"Error adding reaction: {e}")
            return {"success": False, "error": str(e)}

    def update(self, model: BaseNode) -> Dict[str, Any]:
        """Update an existing node in Neo4j."""
        logger.info(f"Updating {model.__class__.__name__} in Neo4j")
        try:
            unique_value = model.get_unique_id()
            unique_field = model.get_unique_field_name()
            logger.debug(f"Model unique_id: {unique_value}")
            logger.debug(f"Model unique_field: {unique_field}")

            if not unique_value:
                logger.error("No unique field found for update")
                return {"success": False, "error": "No unique field found for update"}

            with self.driver.session() as session:
                # Find existing node
                label = model.__class__.__name__
                logger.debug(
                    f"Looking for node with label: {label}, unique_field: {unique_field}, value: {unique_value}"
                )

                cypher = f"MATCH (n:{label}) WHERE n.{unique_field} = $value RETURN n"
                logger.debug(f"Update cypher: {cypher}")

                result = session.run(cypher, {"value": unique_value})
                existing_node = result.single()

                if not existing_node:
                    logger.error(
                        f"Node not found for {label} with {unique_field} = {unique_value}"
                    )
                    return {"success": False, "error": f"Node not found for {label}"}

                logger.debug(f"Found existing node: {existing_node}")

                # Update properties
                properties = self.mapper._create_node_properties(model)
                properties_str = ", ".join([f"n.{k} = ${k}" for k in properties.keys()])
                logger.debug(f"Update properties: {properties}")

                cypher = f"MATCH (n:{label}) WHERE n.{unique_field} = $value SET {properties_str}"
                logger.debug(f"Update SET cypher: {cypher}")

                session.run(cypher, {**properties, "value": unique_value})
                logger.info(f"Successfully updated {model.__class__.__name__}")

                return {"success": True, "node_updated": True}

        except Exception as e:
            logger.error(f"Error updating {model.__class__.__name__} in Neo4j: {e}")
            logger.exception("Full traceback:")
            return {"success": False, "error": str(e)}

    def setup_database(self) -> None:
        """Set up the Neo4j database with constraints and indexes."""
        logger.info("Setting up Neo4j database with constraints and indexes")
        try:
            self.mapper.create_constraints_and_indexes()
            logger.info("Database setup completed successfully")
        except Exception as e:
            logger.error(f"Error setting up database: {e}")
            logger.exception("Full traceback:")


# Example usage
if __name__ == "__main__":
    logger.info("Starting Neo4j integration test")
    # Create a test protein
    import asyncio

    from pyeed.fetch.uniprot import get_proteins_from_uniprot

    from .model import Embedding

    ids = ["P07486", "P69905", "P12345", "P00360"]

    proteins = asyncio.run(get_proteins_from_uniprot(ids))

    # Add the new embedding to the protein BEFORE creating it in Neo4j
    new_embed = Embedding(
        description="esm2_mean_pooled",
        model_name="esm2-t33-650M-UR50S",
        pooling_method="mean",
        vector=[0.1, 0.2, 0.3],
        n_dims=3,
    )
    # Find the P69905 protein specifically
    p69905_protein = next((p for p in proteins if p.accession_id == "P69905"), None)
    if p69905_protein:
        p69905_protein.add_embedding("P69905", new_embed)
    else:
        proteins[0].add_embedding(proteins[0].accession_id, new_embed)

    bd_cred = {
        "uri": "bolt://127.0.0.1:7687",
        "user": "neo4j",
        "password": "12345678",
    }

    writer = PyeedNeo4jWriter(**bd_cred)
    mapper = PyeedNeo4jMapper(writer.driver)
    mapper.create_constraints_and_indexes()

    # Now add proteins using the new add_protein method
    results = []
    for protein in proteins:
        results.append(writer.add_protein(protein))
    print("Main creation result:", results)

    writer.close()
