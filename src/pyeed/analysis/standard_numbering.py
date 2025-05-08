"""Module: standard_numbering
This module provides the StandardNumberingTool class which allows the creation and management
of standard numbering schemes for protein or DNA sequences. It supports aligning sequences
using ClustalOmega for a multiple sequence alignment or using a pairwise alignment method,
and then updates a Neo4j database with the computed numbering positions.

Dependencies:
- loguru: For logging.
- pyeed: For database connection, models, and alignment tools.
- ClustalOmega: For multiple sequence alignment.
"""

from typing import Any, Dict, List, Optional, Tuple

from loguru import logger
from pyeed.analysis.sequence_alignment import PairwiseAligner
from pyeed.dbconnect import DatabaseConnector
from pyeed.model import StandardNumbering
from pyeed.tools.clustalo import ClustalOmega


class StandardNumberingTool:
    """A tool for defining and managing standard numbering schemes for proteins or DNA sequences.

    The standard numbering scheme defines a consistent position numbering system across multiple
    sequences. Each sequence in the database can be linked to a standard numbering definition,
    with the relationship storing an array of positions that map the sequence positions to
    the standard numbering scheme.
    """

    def __init__(self, name: str) -> None:
        """
        Initialize the StandardNumberingTool with a specific scheme name.

        Args:
            name: The name of the standard numbering scheme.

        Attributes:
            positions: A dictionary mapping protein accession ids to lists of numbering positions.
        """
        self.name = name

    def get_node_base_sequence(
        self,
        base_sequence_id: str,
        db: DatabaseConnector,
        node_type: str = "Protein",
        region_ids_neo4j: Optional[list[str]] = None,
    ) -> dict[str, str]:
        """
        Retrieve the base node sequence from the database for a given accession id.

        This method executes a query that returns the node with the provided id. It assumes a
        valid result is returned and constructs a dictionary containing the id and the sequence.

        Args:
            base_sequence_id: The accession id of the base node sequence.
            db: The database connector instance to perform the query.
            region_ids_neo4j: A list of region IDs for the sequence cuting based on region_based_sequence.
        Returns:
            A dictionary with keys 'id' and 'sequence' holding the node type id and its sequence.
        """
        if region_ids_neo4j:
            query = f"""
            MATCH (p:{node_type})-[e:HAS_REGION]->(r:Region)
            WHERE id(r) IN $region_ids_neo4j
            WHERE p.accession_id = '{base_sequence_id}'
            RETURN p.accession_id AS accession_id, e.start AS start, e.end AS end, p.sequence AS sequence
            """
        else:
            query = f"""
            MATCH (p:{node_type})
            WHERE p.accession_id = '{base_sequence_id}'
            RETURN p.accession_id AS accession_id, p.sequence AS sequence
            """
        base_sequence_read = db.execute_read(query)
        # Assume the first returned record is the desired base sequence
        if region_ids_neo4j:
            base_sequence = {
                "id": base_sequence_read[0]["accession_id"],
                "sequence": base_sequence_read[0]["sequence"][
                    base_sequence_read[0]["start"] : base_sequence_read[0]["end"]
                ],
            }
        else:
            base_sequence = {
                "id": base_sequence_read[0]["accession_id"],
                "sequence": base_sequence_read[0]["sequence"],
            }
        return base_sequence

    def save_positions(
        self,
        db: DatabaseConnector,
        positions: dict[str, list[str]],
        node_type: str = "Protein",
        region_ids_neo4j: Optional[list[str]] = None,
    ) -> None:
        """
        Save the calculated numbering positions for each protein into the database.

        For each protein id stored in self.positions, this method creates (or merges) the relationship
        between the protein and its standard numbering node. The positions list for each protein is converted
        to a string (this could be refined later if a different data model is needed).

        Args:
            db: The database connector instance used to execute the write queries.
            positions: A dictionary mapping protein accession ids to lists of numbering positions.
            node_type: The type of node to process. Default is "Protein".
            region_ids_neo4j: A list of region IDs for the sequence cuting based on region_based_sequence.
        """
        for protein_id in positions:
            if region_ids_neo4j:
                query = f"""
                    MATCH (p:{node_type} {{accession_id: '{protein_id}'}})-[e:HAS_REGION]->(r:Region)
                    WHERE id(r) IN $region_ids_neo4j
                    MATCH (s:StandardNumbering {{name: '{self.name}'}})
                    MERGE (r)-[rel:HAS_STANDARD_NUMBERING]->(s)
                    SET rel.positions = {str(positions[protein_id])}
                """
                db.execute_write(
                    query, parameters={"region_ids_neo4j": region_ids_neo4j}
                )
            else:
                query = f"""
                    MATCH (p:{node_type} {{accession_id: '{protein_id}'}})
                    MATCH (s:StandardNumbering {{name: '{self.name}'}})
                    MERGE (p)-[rel:HAS_STANDARD_NUMBERING]->(s)
                    SET rel.positions = {str(positions[protein_id])}
                """
                db.execute_write(query)

    def run_numbering_algorithm_clustalo(
        self, base_sequence_id: str, alignment: Any
    ) -> Dict[str, List[str]]:
        """
        Compute standard numbering positions for each sequence based on a multiple alignment.

        The method uses the first sequence (assumed to be the base) as a reference and assigns
        sequential numbering for every non-gap character. When an insert (a gap in the base) occurs,
        the algorithm assigns a number with an insert notation (e.g. "0.1", "2.1"). Temporarily, gaps
        in target sequences are marked with ".GAP" and are later removed from the final positions.

        Args:
            base_sequence_id: The accession id of the base sequence.
            alignment: The multiple alignment object where each element has:
                       - .seq : the aligned sequence string
                       - .id  : the protein id

        Returns:
            A dictionary mapping each protein's accession id to a list of computed numbering positions as strings.
        """
        logger.info(f"Running numbering algorithm for base sequence {base_sequence_id}")

        positions: Dict[str, List[str]] = {}

        # Convert alignment to a list so that it is subscriptable.
        alignment_list = list(alignment)[0][1]
        # Extract base sequence (first in alignment)
        base_sequence = alignment_list[0].sequence
        # Populate base numbering: count only non-gap characters
        positions[base_sequence_id] = [
            str(i + 1) for i in range(len(base_sequence.replace("-", "")))
        ]

        base_seq_counter = -1
        # Iterate over each aligned position
        for pos in range(len(alignment_list[0].sequence)):
            base_aa = alignment_list[0].sequence[pos]
            # Move counter when a real amino acid is encountered
            if base_aa != "-":
                base_seq_counter += 1

            # Process each target sequence in the alignment (skip base sequence at index 0)
            for j in range(1, len(alignment_list)):
                sequence = alignment_list[j].sequence
                sequence_id = alignment_list[j].id

                # Both sequences have gaps at this position; nothing to number
                if base_aa == "-" and sequence[pos] == "-":
                    continue

                # Handle insert in the base sequence:
                if base_aa == "-" and sequence[pos] != "-":
                    if sequence_id not in positions:
                        positions[sequence_id] = []
                        positions[sequence_id].append("0.1")
                    else:
                        # If the last number is already an insert, increase the insert counter
                        if "." in positions[sequence_id][-1]:
                            insert_number = int(
                                positions[sequence_id][-1].split(".")[1]
                            )
                            base_number = int(positions[sequence_id][-1].split(".")[0])
                            positions[sequence_id].append(
                                f"{base_number}.{insert_number + 1}"
                            )
                        else:
                            positions[sequence_id].append(
                                f"{positions[sequence_id][-1]}.1"
                            )

                # Handle deletion in the target sequence (i.e., gap in target):
                if base_aa != "-" and sequence[pos] == "-":
                    if sequence_id not in positions:
                        # Mark gap with a .GAP suffix based on corresponding base position
                        positions[sequence_id] = [
                            positions[base_sequence_id][base_seq_counter] + ".GAP"
                        ]
                    positions[sequence_id].append(
                        f"{positions[sequence_id][-1].split('.')[0]}.GAP"
                    )

                # Both positions are amino acids; advance numbering normally.
                if base_aa != "-" and sequence[pos] != "-":
                    if sequence_id not in positions:
                        positions[sequence_id] = []
                        positions[sequence_id].append(
                            positions[base_sequence_id][base_seq_counter]
                        )
                    else:
                        # If previous number was an insert, increment the numerical part
                        if (
                            "." in positions[sequence_id][-1]
                            and "GAP" not in positions[sequence_id][-1]
                        ):
                            base_number = int(positions[sequence_id][-1].split(".")[0])
                            positions[sequence_id].append(f"{base_number + 1}")
                        else:
                            positions[sequence_id].append(
                                f"{int(positions[base_sequence_id][base_seq_counter])}"
                            )

        # Remove any gap placeholders before returning the positions.
        for protein_id in positions:
            positions[protein_id] = [
                pos
                for pos in positions[protein_id]
                if pos is not None and "GAP" not in pos
            ]

        return positions

    def run_numbering_algorithm_pairwise(
        self, base_sequence_id: str, alignment: List[Tuple[str, str, str]]
    ) -> Dict[str, List[str]]:
        """
        Compute numbering positions using pairwise alignments relative to a base sequence.

        The alignment is expected to be a list where each element is a tuple containing:
            (base_sequence_aligned, target_sequence_aligned, target_protein_id)
        The function assigns standard numbering by considering inserts and deletions as in the
        multiple alignment approach, and finally removes any gap placeholders.

        Args:
            base_sequence_id: The accession id of the base sequence.
            alignment: A list of tuples with (base_sequence, target_sequence, target_id).

        Returns:
            A dictionary mapping protein accession ids to their list of numbering positions.
        """
        positions: Dict[str, List[str]] = {}
        # Initialize base positions by counting non-gap characters in the base aligned sequence.
        positions[base_sequence_id] = [
            str(i + 1) for i in range(len(alignment[0][0].replace("-", "")))
        ]

        # Process each pairwise alignment tuple.
        for pair in alignment:
            base_sequence = pair[0]
            target_sequence = pair[1]
            target_sequence_id = pair[2]

            base_seq_counter = -1

            for pos in range(len(base_sequence)):
                # Increment counter for non-gap in base sequence
                if base_sequence[pos] != "-":
                    base_seq_counter += 1

                # Skip positions where both sequences have a gap.
                if base_sequence[pos] == "-" and target_sequence[pos] == "-":
                    continue

                # Handle insert in the base sequence.
                if base_sequence[pos] == "-" and target_sequence[pos] != "-":
                    if target_sequence_id not in positions:
                        positions[target_sequence_id] = []
                        positions[target_sequence_id].append("0.1")
                    else:
                        if "." in positions[target_sequence_id][-1]:
                            insert_number = int(
                                positions[target_sequence_id][-1].split(".")[1]
                            )
                            base_number = int(
                                positions[target_sequence_id][-1].split(".")[0]
                            )
                            positions[target_sequence_id].append(
                                f"{base_number}.{insert_number + 1}"
                            )
                        else:
                            positions[target_sequence_id].append(
                                f"{positions[target_sequence_id][-1]}.1"
                            )

                # Handle deletion in the target sequence.
                if base_sequence[pos] != "-" and target_sequence[pos] == "-":
                    if target_sequence_id not in positions:
                        positions[target_sequence_id] = []
                        positions[target_sequence_id].append(
                            f"{positions[base_sequence_id][pos]}.GAP"
                        )
                    positions[target_sequence_id].append(
                        f"{positions[target_sequence_id][-1].split('.')[0]}.GAP"
                    )

                # Both positions are amino acids; update numbering regularly.
                if base_sequence[pos] != "-" and target_sequence[pos] != "-":
                    if target_sequence_id not in positions:
                        positions[target_sequence_id] = []
                        positions[target_sequence_id].append(
                            positions[base_sequence_id][base_seq_counter]
                        )
                    else:
                        if (
                            "." in positions[target_sequence_id][-1]
                            and "GAP" not in positions[target_sequence_id][-1]
                        ):
                            base_number = int(
                                positions[target_sequence_id][-1].split(".")[0]
                            )
                            positions[target_sequence_id].append(f"{base_number + 1}")
                        else:
                            positions[target_sequence_id].append(
                                f"{int(positions[base_sequence_id][base_seq_counter])}"
                            )

        # Remove gap placeholders before returning
        for protein_id in positions:
            positions[protein_id] = [
                pos
                for pos in positions[protein_id]
                if pos is not None and "GAP" not in pos
            ]

        return positions

    def apply_standard_numbering_pairwise(
        self,
        base_sequence_id: str,
        db: DatabaseConnector,
        list_of_seq_ids: Optional[List[str]] = None,
        return_positions: bool = False,
        node_type: str = "Protein",
        region_ids_neo4j: Optional[list[str]] = None,
    ) -> Optional[Dict[str, List[str]]]:
        """
        Apply standard numbering via pairwise alignment using a base sequence.

        If a specific list of sequence ids is not provided, this function fetches all sequence ids
        from the database. It then creates pairs with the base sequence as the reference and performs
        a pairwise alignment using the PairwiseAligner. The computed numbering positions are stored,
        a StandardNumbering node is created (or retrieved), and the positions are saved to the database.

        Args:
            base_sequence_id: The accession id of the base sequence.
            db: The DatabaseConnector instance used for communication with the database.
            list_of_seq_ids: An optional list of node type ids to process. If None, all node type ids are used.
            return_positions: If True, the method returns the computed positions dictionary after processing.
            node_type: The type of node to process. Default is "Protein".
            region_ids_neo4j: A list of region IDs for the sequence cuting based on region_based_sequence.
        Raises:
            ValueError: If the pairwise alignment fails and returns no results.
        """
        if list_of_seq_ids is None:
            query = f"""
            MATCH (p:{node_type})
            WHERE p.accession_id IS NOT NULL
            RETURN p.accession_id AS accession_id
            """
            results = db.execute_read(query)
            if results is None:
                raise ValueError("No results returned from the query")
            for row in results:
                if row is None:
                    raise ValueError("No results returned from the query")
                if not row.get("accession_id"):
                    raise ValueError("Row missing required accession_id field")
            list_of_seq_ids = [row["accession_id"] for row in results]

        # Remove the base sequence id from the list if present.
        while base_sequence_id in list_of_seq_ids:
            list_of_seq_ids.remove(base_sequence_id)

        # Generate pairs with the base sequence as the first element.
        pairs = []
        for node_id in list_of_seq_ids:
            pairs.append((base_sequence_id, node_id))

        # check if the pairs are already existing with the same name under the same standard numbering node
        if node_type == "DNA" and region_ids_neo4j is not None:
            query = """
            MATCH (s:StandardNumbering {name: $name})
            MATCH (d:DNA)-[e:HAS_REGION]-(r:Region)-[:HAS_STANDARD_NUMBERING]-(s)
            WHERE id(r) IN $region_ids_neo4j
            AND d.accession_id IN $list_of_seq_ids
            RETURN d.accession_id AS accession_id
            """

            results = db.execute_read(
                query,
                parameters={
                    "list_of_seq_ids": list_of_seq_ids,
                    "name": self.name,
                    "region_ids_neo4j": region_ids_neo4j,
                },
            )
        else:
            query = f"""
            MATCH (s:StandardNumbering {{name: $name}})
            MATCH (p:{node_type})-[rel:HAS_STANDARD_NUMBERING]->(s)
            WHERE p.accession_id IN $list_of_seq_ids
            RETURN p.accession_id AS accession_id
            """
            results = db.execute_read(
                query,
                parameters={"list_of_seq_ids": list_of_seq_ids, "name": self.name},
            )

        if results is not None:
            for row in results:
                if row is not None:
                    if row.get("accession_id"):
                        pairs.remove((base_sequence_id, row["accession_id"]))
                        logger.info(
                            f"Pair {base_sequence_id} and {row['accession_id']} already exists under the same standard numbering node"
                        )

        # remove double pairs in the list of pairs
        pairs = list(set(pairs))
        logger.info(f"Pairs: {pairs}")

        # Run the pairwise alignment using the PairwiseAligner.
        pairwise_aligner = PairwiseAligner(node_type=node_type)
        input = (list_of_seq_ids or []) + [base_sequence_id]
        if not input:
            raise ValueError("No input sequences provided")

        logger.info(f"Input: {input}")

        results_pairwise = pairwise_aligner.align_multipairwise(
            ids=input,  # Combine ids for alignment
            db=db,
            pairs=pairs,  # List of sequence pairs to be aligned
            node_type=node_type,
            region_ids_neo4j=region_ids_neo4j,
        )

        # logger.info(f"Pairwise alignment results: {results_pairwise}")

        if results_pairwise is None:
            raise ValueError("Pairwise alignment failed - no results returned")

        # Convert results from dicts to a list of tuples with desired order.
        converted_alignment: List[Tuple[str, str, str]] = [
            (
                str(result["query_aligned"]),
                str(result["target_aligned"]),
                str(result["target_id"]),
            )
            for result in results_pairwise
        ]

        if len(converted_alignment) == 0:
            logger.info(f"No alignment found for {base_sequence_id}")
            return None

        logger.info(f"Converted alignment: {len(converted_alignment)}")

        # Compute positions using the pairwise numbering algorithm.
        positions = self.run_numbering_algorithm_pairwise(
            base_sequence_id, converted_alignment
        )

        # Ensure the standard numbering node exists in the database.
        StandardNumbering.get_or_save(
            name=self.name,
            definition=f"Pairwise based on base sequence {base_sequence_id}",
        )

        # Update the database with the calculated positions.
        self.save_positions(db, positions, node_type, region_ids_neo4j)

        if return_positions:
            return positions
        return None

    def apply_standard_numbering(
        self,
        base_sequence_id: str,
        db: DatabaseConnector,
        list_of_seq_ids: Optional[List[str]] = None,
        node_type: str = "Protein",
        region_ids_neo4j: Optional[list[str]] = None,
    ) -> None:
        """
        Apply a standard numbering scheme to a collection of nodes using multiple sequence alignment.

        This method first retrieves all node sequences from the database (or a subset if list_of_seq_ids is provided).
        It then uses ClustalOmega to perform a multiple sequence alignment, computes the numbering positions via
        run_numbering_algorithm, creates (or retrieves) a StandardNumbering node, and saves the positions back into the database.

        Args:
            base_sequence_id: The accession id of the base sequence to which others are aligned.
            db: DatabaseConnector instance used for executing queries.
            list_of_seq_ids: An optional list of specific node type ids to process. If None, all node type ids are used.
            node_type: The type of node to process. Default is "Protein".
            region_ids_neo4j: A list of region IDs for the sequence cuting based on region_based_sequence.
        """

        if list_of_seq_ids is None:
            query = f"""
            MATCH (p:{node_type}) 
            WHERE p.sequence IS NOT NULL
            RETURN p.accession_id AS accession_id
            """
            results = db.execute_read(query)
            if results is None:
                raise ValueError("No results returned from the query")
            list_of_seq_ids = [row["accession_id"] for row in results]

        # Retrieve all nodes from the database. With both id and sequence.
        query = f"""
        MATCH (p:{node_type})
        WHERE p.sequence IS NOT NULL
        AND p.accession_id IN $list_of_seq_ids
        RETURN p.accession_id AS accession_id, p.sequence AS sequence
        """
        # Execute the query and build the nodes dictionary
        nodes_read: List[Dict[str, Any]]
        query_result = db.execute_read(
            query, parameters={"list_of_seq_ids": list_of_seq_ids}
        )
        if query_result is None:
            nodes_read = []
        else:
            nodes_read = query_result

        if node_type == "DNA" and region_ids_neo4j is not None:
            # then the sequence is a region based sequence.
            # get the region objects for each of the nodes as well
            query = f"""
            MATCH (p:{node_type})-[e:HAS_REGION]->(r:Region)
            WHERE id(r) IN $region_ids_neo4j
            WHERE p.accession_id IN $list_of_seq_ids
            RETURN p.accession_id AS accession_id, e.start AS start, e.end AS end, p.sequence AS sequence
            """
            region_read = db.execute_read(
                query,
                parameters={
                    "list_of_seq_ids": list_of_seq_ids,
                    "region_ids_neo4j": region_ids_neo4j,
                },
            )
            nodes_dict = {
                node["accession_id"]: node["sequence"][node["start"] : node["end"]]
                for node in region_read
            }

        else:
            nodes_dict = {node["accession_id"]: node["sequence"] for node in nodes_read}

        logger.info(f"Using {len(nodes_dict)} sequences for standard numbering")

        # Obtain the base sequence details from the database.
        base_sequence = self.get_node_base_sequence(
            base_sequence_id, db, node_type, region_ids_neo4j
        )

        # Remove the base sequence from the nodes list to prevent duplicate alignment.
        if base_sequence_id in nodes_dict:
            nodes_dict.pop(base_sequence_id)

        # Create a dictionary for ClustalOmega that includes both the base and target sequences.
        sequences_dict = {base_sequence["id"]: base_sequence["sequence"]}
        for key in nodes_dict:
            sequences_dict[key] = nodes_dict[key]

        # Run the multiple sequence alignment using ClustalOmega.
        clustalO = ClustalOmega()
        alignment = clustalO.align(
            sequences_dict
        )  # Passing a dict of sequences to ClustalOmega.

        logger.info(f"Alignment received from ClustalOmega:\n{alignment}")
        logger.info(f"Alignment length: {len(list(alignment)[0][1][0].sequence)}")

        # Compute standard numbering positions using the computed alignment.
        positions = self.run_numbering_algorithm_clustalo(base_sequence_id, alignment)
        logger.info(f"Positions computed: {positions}")

        # Create (or get) the StandardNumbering node in the database.
        StandardNumbering.get_or_save(
            name=f"{self.name}",
            definition=f"ClustalO based on base sequence {base_sequence_id}",
        )

        # Update the database with the relationships between nodes and standard numbering.
        self.save_positions(db, positions, node_type, region_ids_neo4j)
