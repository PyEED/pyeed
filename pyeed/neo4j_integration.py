"""Generic Neo4j integration for Pyeed models using NodeHint and EdgeHint metadata."""

import logging
from enum import Enum
from typing import Any, Dict, List, Optional, Set, TypeVar

from neo4j import GraphDatabase
from pydantic import BaseModel

from .model import BaseNode, EdgeHint, NodeHint

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


def extract_edge_hints(model: type[BaseModel]) -> Dict[str, EdgeHint]:
    """Extract EdgeHint metadata from model fields."""
    hints: Dict[str, EdgeHint] = {}
    logger.debug(f"Extracting EdgeHints from model: {model.__name__}")

    for field_name, field_info in model.model_fields.items():
        for meta in getattr(field_info, "metadata", ()):
            if isinstance(meta, EdgeHint):
                hints[field_name] = meta
                logger.debug(f"Found EdgeHint for field '{field_name}': {meta}")
                break

    logger.debug(f"Extracted {len(hints)} EdgeHints: {hints}")
    return hints


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

    def create_relationships(
        self, source_model: BaseNode, target_models: List[BaseNode], edge_hint: EdgeHint
    ) -> List[Neo4jRelationship]:
        """Create relationships based on EdgeHint metadata."""
        logger.debug(
            f"Creating relationships from {source_model.__class__.__name__} to {len(target_models)} targets"
        )
        logger.debug(f"EdgeHint: {edge_hint}")

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

            # Create relationship properties
            rel_properties = {
                "source_key": getattr(source_model, edge_hint.outgoing_attr_name, None),
                "target_key": getattr(target_model, edge_hint.outgoing_attr_name, None),
            }
            logger.debug(f"Relationship properties: {rel_properties}")

            rel = Neo4jRelationship(
                start_node_id=source_id,
                end_node_id=target_id,
                rel_type=edge_hint.name.upper(),
                properties=rel_properties,
            )
            logger.debug(f"Created relationship: {rel}")
            relationships.append(rel)

        logger.debug(f"Created {len(relationships)} relationships")
        return relationships

    def map_model_to_neo4j(self, model: BaseNode) -> Dict[str, Any]:
        """Map any BaseNode model to Neo4j nodes and relationships."""
        logger.debug(f"Mapping model to Neo4j: {model.__class__.__name__}")
        result: Dict[str, Any] = {"nodes": [], "relationships": []}

        # Create main model node
        logger.debug("Creating main model node...")
        main_node = self.create_node(model)
        result["nodes"].append(main_node)
        logger.debug(f"Added main node: {main_node}")

        # Check ALL fields for nested BaseNode objects, not just those with EdgeHints
        logger.debug("Scanning all fields for nested BaseNode objects...")

        for field_name, field_info in type(model).model_fields.items():
            field_value = getattr(model, field_name, None)
            logger.debug(f"Checking field '{field_name}': {type(field_value)}")

            if field_value is None:
                logger.debug(f"Field {field_name} is None, skipping")
                continue

            # Check for lists of BaseNode objects
            if isinstance(field_value, list) and field_value:
                logger.debug(
                    f"Found list field: {field_name} with {len(field_value)} items"
                )

                if all(isinstance(item, BaseNode) for item in field_value):
                    logger.debug(f"All items in {field_name} are BaseNode instances")

                    # Create nodes for each related model
                    for i, related_model in enumerate(field_value):
                        logger.debug(
                            f"Creating node {i+1} for related model: {related_model.__class__.__name__}"
                        )
                        related_node = self.create_node(related_model)
                        result["nodes"].append(related_node)
                        logger.debug(f"Added related node: {related_node}")

                    # Check if there's an EdgeHint for this field to create relationships
                    edge_hints = extract_edge_hints(type(model))
                    if field_name in edge_hints:
                        edge_hint = edge_hints[field_name]
                        logger.debug(
                            f"Creating relationships for list field: {field_name} using EdgeHint"
                        )
                        relationships = self.create_relationships(
                            model, field_value, edge_hint
                        )
                        result["relationships"].extend(relationships)
                        logger.debug(f"Added {len(relationships)} relationships")
                    else:
                        # Use EdgeHints from the satellite models (reverse direction)
                        logger.debug(
                            "No EdgeHint in parent, checking satellite models for reverse relationships"
                        )
                        for related_model in field_value:
                            satellite_edge_hints = extract_edge_hints(
                                type(related_model)
                            )
                            for (
                                sat_field_name,
                                sat_edge_hint,
                            ) in satellite_edge_hints.items():
                                if (
                                    sat_edge_hint.outgoing_class_name
                                    == model.__class__.__name__
                                ):
                                    logger.debug(
                                        f"Found reverse EdgeHint: {sat_edge_hint}"
                                    )
                                    # Create reverse relationship: parent -> satellite
                                    # The EdgeHint in the satellite points TO the parent
                                    # So we create: parent -> satellite
                                    rel_properties = {
                                        "source_key": getattr(
                                            model,
                                            sat_edge_hint.outgoing_attr_name,
                                            None,
                                        ),
                                        "target_key": getattr(
                                            related_model,
                                            sat_edge_hint.outgoing_attr_name,
                                            None,
                                        ),
                                    }
                                    rel = Neo4jRelationship(
                                        start_node_id=model.get_unique_id(),
                                        end_node_id=related_model.get_unique_id(),
                                        rel_type=sat_edge_hint.name.upper(),
                                        properties=rel_properties,
                                    )
                                    result["relationships"].append(rel)
                                    logger.debug(f"Added reverse relationship: {rel}")

            # Check for dicts of BaseNode objects
            elif isinstance(field_value, dict) and field_value:
                logger.debug(
                    f"Found dict field: {field_name} with {len(field_value)} items"
                )

                if all(isinstance(item, BaseNode) for item in field_value.values()):
                    logger.debug(f"All values in {field_name} are BaseNode instances")

                    # Create nodes for each related model
                    for key, related_model in field_value.items():
                        logger.debug(
                            f"Creating node for key '{key}': {related_model.__class__.__name__}"
                        )
                        related_node = self.create_node(related_model)
                        result["nodes"].append(related_node)
                        logger.debug(f"Added related node: {related_node}")

                    # Check if there's an EdgeHint for this field to create relationships
                    edge_hints = extract_edge_hints(type(model))
                    if field_name in edge_hints:
                        edge_hint = edge_hints[field_name]
                        logger.debug(
                            f"Creating relationships for dict field: {field_name} using EdgeHint"
                        )
                        relationships = self.create_relationships(
                            model, list(field_value.values()), edge_hint
                        )
                        result["relationships"].extend(relationships)
                        logger.debug(f"Added {len(relationships)} relationships")
                    else:
                        # Use EdgeHints from the satellite models (reverse direction)
                        logger.debug(
                            "No EdgeHint in parent, checking satellite models for reverse relationships"
                        )
                        for related_model in field_value.values():
                            satellite_edge_hints = extract_edge_hints(
                                type(related_model)
                            )
                            for (
                                sat_field_name,
                                sat_edge_hint,
                            ) in satellite_edge_hints.items():
                                if (
                                    sat_edge_hint.outgoing_class_name
                                    == model.__class__.__name__
                                ):
                                    logger.debug(
                                        f"Found reverse EdgeHint: {sat_edge_hint}"
                                    )
                                    # Create reverse relationship: parent -> satellite
                                    # The EdgeHint in the satellite points TO the parent
                                    # So we create: parent -> satellite
                                    rel_properties = {
                                        "source_key": getattr(
                                            model,
                                            sat_edge_hint.outgoing_attr_name,
                                            None,
                                        ),
                                        "target_key": getattr(
                                            related_model,
                                            sat_edge_hint.outgoing_attr_name,
                                            None,
                                        ),
                                    }
                                    rel = Neo4jRelationship(
                                        start_node_id=model.get_unique_id(),
                                        end_node_id=related_model.get_unique_id(),
                                        rel_type=sat_edge_hint.name.upper(),
                                        properties=rel_properties,
                                    )
                                    result["relationships"].append(rel)
                                    logger.debug(f"Added reverse relationship: {rel}")

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

    def create(self, model: BaseNode) -> Dict[str, Any]:
        """Create a new node in Neo4j for any BaseNode model with automatic relationship handling."""
        logger.info(f"Creating {model.__class__.__name__} in Neo4j")
        try:
            with self.driver.session() as session:
                # First, create the main node
                logger.debug(f"Creating main node for {model.__class__.__name__}")
                neo4j_node = self.mapper.create_node(model)

                # Create/merge the node in Neo4j using unique constraint
                labels_str = ":".join(neo4j_node.labels)
                unique_field = model.get_unique_field_name()
                unique_value = model.get_unique_id()

                # Use MERGE on unique field to avoid duplicates
                cypher = f"""
                MERGE (n:{labels_str} {{{unique_field}: $unique_value}})
                SET n += $properties
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
                    logger.error(f"Failed to create node: {neo4j_node}")
                    return {"success": False, "error": "Failed to create node"}

                neo4j_id = str(created_node["n"].element_id)
                logger.info(
                    f"Created node {neo4j_node.labels} with Neo4j ID {neo4j_id}"
                )

                # Check if this model has EdgeHints that point to existing nodes
                edge_hints = extract_edge_hints(type(model))
                relationships_created = 0

                if edge_hints:
                    logger.debug(
                        f"Found {len(edge_hints)} EdgeHints, checking for existing target nodes"
                    )

                    for field_name, edge_hint in edge_hints.items():
                        logger.debug(f"Processing EdgeHint for field: {field_name}")

                        # Get the value that should link to the target
                        source_value = getattr(
                            model, edge_hint.outgoing_attr_name, None
                        )
                        if not source_value:
                            logger.debug(
                                f"No value for outgoing_attr_name {edge_hint.outgoing_attr_name}, skipping"
                            )
                            continue

                        logger.debug(
                            f"Looking for existing {edge_hint.outgoing_class_name} with {edge_hint.outgoing_attr_name} = {source_value}"
                        )

                        # Check if target node exists
                        target_cypher = f"MATCH (target:{edge_hint.outgoing_class_name} {{{edge_hint.outgoing_attr_name}: $target_value}}) RETURN target"
                        target_result = session.run(
                            target_cypher, {"target_value": source_value}
                        )
                        target_node = target_result.single()

                        if target_node:
                            logger.debug(
                                "Found existing target node, creating relationship"
                            )

                            # Create relationship TO the target class using elementId()
                            # The EdgeHint is on the field that connects TO the target
                            rel_cypher = f"""
                            MATCH (target:{edge_hint.outgoing_class_name} {{{edge_hint.outgoing_attr_name}: $target_value}})
                            MATCH (source:{model.__class__.__name__}) WHERE elementId(source) = $source_element_id
                            MERGE (target)-[r:{edge_hint.name}]->(source)
                            SET r += {{source_key: $target_value, target_key: $source_unique_id}}
                            RETURN r
                            """

                            logger.debug(f"Relationship cypher: {rel_cypher}")
                            rel_result = session.run(
                                rel_cypher,
                                {
                                    "target_value": source_value,
                                    "source_element_id": neo4j_id,  # Use the elementId from node creation
                                    "source_unique_id": model.get_unique_id(),
                                },
                            )

                            if rel_result.single():
                                relationships_created += 1
                                logger.info(
                                    f"Created relationship {edge_hint.name} from {model.__class__.__name__} to {edge_hint.outgoing_class_name}"
                                )
                            else:
                                logger.warning(
                                    f"Failed to create relationship for field {field_name}"
                                )
                        else:
                            logger.debug(
                                f"No existing {edge_hint.outgoing_class_name} found with {edge_hint.outgoing_attr_name} = {source_value}"
                            )

                # Also handle any nested models if this is a complex model with lists/dicts of BaseNode objects
                nested_nodes_created = 0
                nested_relationships_created = 0

                # Map the model to check for nested structures
                neo4j_structure = self.mapper.map_model_to_neo4j(model)
                if len(neo4j_structure["nodes"]) > 1:  # More than just the main node
                    logger.debug(
                        f"Found {len(neo4j_structure['nodes']) - 1} nested nodes to create"
                    )

                    # Create nested nodes (skip the first one which is the main node we already created)
                    node_id_mapping = {neo4j_node.unique_id: neo4j_id}

                    for nested_node in neo4j_structure["nodes"][1:]:
                        nested_labels_str = ":".join(nested_node.labels)

                        # Find unique field for nested node (we need to determine its class)
                        # Extract class name from labels (assume single label per node)
                        nested_class_name = list(nested_node.labels)[0]
                        nested_unique_field = "_id"  # Default fallback
                        nested_unique_value = nested_node.unique_id

                        # Try to get proper unique field from class
                        try:
                            import sys

                            from .model import BaseNode

                            base_module = sys.modules.get("pyeed.model")
                            if base_module:
                                nested_model_class = getattr(
                                    base_module, nested_class_name, None
                                )
                                if nested_model_class and issubclass(
                                    nested_model_class, BaseNode
                                ):
                                    nested_unique_field = nested_model_class.get_unique_field_name_for_class()
                        except Exception as e:
                            logger.debug(
                                f"Could not determine unique field for {nested_class_name}, using _id: {e}"
                            )

                        nested_cypher = f"""
                        MERGE (n:{nested_labels_str} {{{nested_unique_field}: $unique_value}})
                        SET n += $properties
                        RETURN n
                        """

                        nested_result = session.run(
                            nested_cypher,
                            {
                                "unique_value": nested_unique_value,
                                "properties": nested_node.properties,
                            },
                        )
                        nested_created = nested_result.single()

                        if nested_created:
                            nested_neo4j_id = str(nested_created["n"].element_id)
                            node_id_mapping[nested_node.unique_id] = nested_neo4j_id
                            nested_nodes_created += 1
                            logger.info(
                                f"Created nested node {nested_node.labels} with Neo4j ID {nested_neo4j_id}"
                            )

                    # Create relationships for nested structures
                    for rel in neo4j_structure["relationships"]:
                        start_neo4j_id = node_id_mapping.get(rel.start_node_id)
                        end_neo4j_id = node_id_mapping.get(rel.end_node_id)

                        if start_neo4j_id and end_neo4j_id:
                            rel_cypher = (
                                f"MATCH (a), (b) "
                                f"WHERE elementId(a) = $start_id AND elementId(b) = $end_id "
                                f"MERGE (a)-[r:{rel.rel_type}]->(b) "
                                f"SET r += $properties "
                                f"RETURN r"
                            )

                            session.run(
                                rel_cypher,
                                {
                                    "start_id": start_neo4j_id,
                                    "end_id": end_neo4j_id,
                                    "properties": rel.properties,
                                },
                            )
                            nested_relationships_created += 1
                            logger.info(f"Created nested relationship {rel.rel_type}")

                total_nodes = 1 + nested_nodes_created
                total_relationships = (
                    relationships_created + nested_relationships_created
                )

                logger.info(
                    f"Successfully created {total_nodes} nodes and {total_relationships} relationships"
                )
                return {
                    "success": True,
                    "nodes_created": total_nodes,
                    "relationships_created": total_relationships,
                }

        except Exception as e:
            logger.error(f"Error creating {model.__class__.__name__} in Neo4j: {e}")
            logger.exception("Full traceback:")
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
    from .model import AnnotationType, Embedding, Organism, Protein, SequenceAnnotation

    # Create a test protein

    prot = Protein(
        accession_id="P01234",
        sequence="MALWMRLLPLLALLALWGPDPAAA",
        name="Test Protein",
        seq_length=24,
        mol_weight=1000,
        ec_number="1.2.3.4",
        nucleotide_id="N01234",
        nucleotide_start=1,
        nucleotide_end=25,
        locus_tag="L01234",
        custom={
            "mentioned_in": [
                "doi:10.1016/j.xinn.2025.100344",
                "doi:10.1016/j.xinn.2025.100345",
            ],
            "already_characterized": False,
            "test": 123,
        },
    )

    prot.add_annotation(
        SequenceAnnotation(
            accession_id="P01234",
            annotation_type=AnnotationType.ACTIVE_SITE,
            positions=[1, 2, 3],
            custom={"my_custom_evidence": "literature", "validated": True},
        )
    )

    prot.add_embedding(
        Embedding(
            accession_id="P01234",
            description="esm2_mean_pooled",
            model_name="esm2-t33-650M-UR50S",
            pooling_method="mean",
            vector=[0.1, 0.2, 0.3],
            n_dims=3,
        ),
    )

    prot.add_organism(
        Organism(
            accession_id="P01234",
            tax_id=123456,
            name="Test Organism",
        )
    )

    # Solution 1: Add new embedding to the existing protein first
    new_embed = Embedding(
        accession_id="P01234",
        description="esm2_mean_pooled_last_hidden_state",
        model_name="esm2-t33-650M-UR50S",
        pooling_method="mean",
        vector=[1.4, 1.5, 1.6],
        n_dims=3,
    )

    # Add the new embedding to the protein BEFORE creating it in Neo4j
    prot.add_embedding(new_embed)

    bd_cred = {
        "uri": "bolt://127.0.0.1:7687",
        "user": "neo4j",
        "password": "12345678",
    }

    writer = PyeedNeo4jWriter(**bd_cred)
    mapper = PyeedNeo4jMapper(writer.driver)
    mapper.create_constraints_and_indexes()

    # Now create everything at once - this will create the relationships automatically
    result = writer.create(prot)
    print("Main creation result:", result)

    # Example of adding a new embedding to an EXISTING protein (simplified!)
    print("\n--- Adding to existing protein ---")
    third_embed = Embedding(
        accession_id="P01234",  # This EdgeHint will automatically link to existing protein
        description="anotherone",
        model_name="UR50S",
        pooling_method="cls",
        vector=[2.1, 2.2, 2.3],
        n_dims=3,
    )

    # Just create the embedding - it will automatically connect to the protein!
    add_result = writer.create(third_embed)
    print("Add embedding result:", add_result)

    writer.close()
