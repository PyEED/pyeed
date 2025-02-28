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
        self.positions: dict[str, list[str]] = {}

    def get_protein_base_sequence(
        self, base_sequence_id: str, db: DatabaseConnector
    ) -> dict[str, str]:
        """
        Retrieve the base protein sequence from the database for a given accession id.

        This method executes a query that returns the protein with the provided id. It assumes a
        valid result is returned and constructs a dictionary containing the id and the sequence.

        Args:
            base_sequence_id: The accession id of the base protein sequence.
            db: The database connector instance to perform the query.

        Returns:
            A dictionary with keys 'id' and 'sequence' holding the protein accession id and its sequence.
        """
        query = f"""
        MATCH (p:Protein)
        WHERE p.accession_id = '{base_sequence_id}'
        RETURN p.accession_id AS accession_id, p.sequence AS sequence
        """
        base_sequence_read = db.execute_read(query)
        # Assume the first returned record is the desired base sequence
        base_sequence = {
            "id": base_sequence_read[0]["accession_id"],
            "sequence": base_sequence_read[0]["sequence"],
        }
        return base_sequence

    def save_positions(self, db: DatabaseConnector) -> None:
        """
        Save the calculated numbering positions for each protein into the database.

        For each protein id stored in self.positions, this method creates (or merges) the relationship
        between the protein and its standard numbering node. The positions list for each protein is converted
        to a string (this could be refined later if a different data model is needed).

        Args:
            db: The database connector instance used to execute the write queries.
        """
        for protein_id in self.positions:
            query = f"""
                MATCH (p:Protein {{accession_id: '{protein_id}'}})
                MATCH (s:StandardNumbering {{name: '{self.name}'}})
                MERGE (p)-[r:HAS_STANDARD_NUMBERING]->(s)
                SET r.positions = {str(self.positions[protein_id])}
            """
            # Execute the write query to update the standard numbering relationship.
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
            list_of_seq_ids: An optional list of protein ids to process. If None, all proteins are used.
            return_positions: If True, the method returns the computed positions dictionary after processing.

        Raises:
            ValueError: If the pairwise alignment fails and returns no results.
        """
        if list_of_seq_ids is None:
            query = """
            MATCH (p:Protein)
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
        list_of_seq_ids.remove(base_sequence_id)

        # Generate pairs with the base sequence as the first element.
        pairs = []
        for protein_id in list_of_seq_ids:
            pairs.append((base_sequence_id, protein_id))

        # Run the pairwise alignment using the PairwiseAligner.
        pairwise_aligner = PairwiseAligner()

        input = (list_of_seq_ids or []) + [base_sequence_id]
        if not input:
            raise ValueError("No input sequences provided")

        results_pairwise = pairwise_aligner.align_multipairwise(
            ids=input,  # Combine ids for alignment
            db=db,
            pairs=pairs,  # List of sequence pairs to be aligned
        )

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

        # Compute positions using the pairwise numbering algorithm.
        self.positions = self.run_numbering_algorithm_pairwise(
            base_sequence_id, converted_alignment
        )

        # Ensure the standard numbering node exists in the database.
        StandardNumbering.get_or_save(
            name=self.name,
            definition=f"Pairwise based on base sequence {base_sequence_id}",
        )

        # Update the database with the calculated positions.
        self.save_positions(db)

        if return_positions:
            return self.positions
        return None

    def apply_standard_numbering(
        self,
        base_sequence_id: str,
        db: DatabaseConnector,
        list_of_seq_ids: Optional[List[str]] = None,
    ) -> None:
        """
        Apply a standard numbering scheme to a collection of proteins using multiple sequence alignment.

        This method first retrieves all protein sequences from the database (or a subset if list_of_seq_ids is provided).
        It then uses ClustalOmega to perform a multiple sequence alignment, computes the numbering positions via
        run_numbering_algorithm, creates (or retrieves) a StandardNumbering node, and saves the positions back into the database.

        Args:
            base_sequence_id: The accession id of the base sequence to which others are aligned.
            db: DatabaseConnector instance used for executing queries.
            list_of_seq_ids: An optional list of specific protein ids to process. If None, all proteins are used.
        """

        if list_of_seq_ids is None:
            query = """
            MATCH (p:Protein) 
            WHERE p.sequence IS NOT NULL
            RETURN p.accession_id AS accession_id
            """
            results = db.execute_read(query)
            if results is None:
                raise ValueError("No results returned from the query")
            list_of_seq_ids = [row["accession_id"] for row in results]

        # Retrieve all proteins from the database. With both id and sequence.
        query = """
        MATCH (p:Protein)
        WHERE p.sequence IS NOT NULL
        AND p.accession_id IN $list_of_seq_ids
        RETURN p.accession_id AS accession_id, p.sequence AS sequence
        """
        # Execute the query and build the proteins dictionary
        proteins_read: List[Dict[str, Any]]
        query_result = db.execute_read(
            query, parameters={"list_of_seq_ids": list_of_seq_ids}
        )
        if query_result is None:
            proteins_read = []
        else:
            proteins_read = query_result
        proteins_dict = {
            protein["accession_id"]: protein["sequence"] for protein in proteins_read
        }

        logger.info(f"Using {len(proteins_dict)} sequences for standard numbering")

        # Obtain the base sequence details from the database.
        base_sequence = self.get_protein_base_sequence(base_sequence_id, db)

        # Remove the base sequence from the proteins list to prevent duplicate alignment.
        proteins_dict.pop(base_sequence_id)

        # Create a dictionary for ClustalOmega that includes both the base and target sequences.
        sequences_dict = {base_sequence["id"]: base_sequence["sequence"]}
        for key in proteins_dict:
            sequences_dict[key] = proteins_dict[key]

        # Run the multiple sequence alignment using ClustalOmega.
        clustalO = ClustalOmega()
        alignment = clustalO.align(
            sequences_dict
        )  # Passing a dict of sequences to ClustalOmega.

        logger.info(f"Alignment received from ClustalOmega:\n{alignment}")

        # Compute standard numbering positions using the computed alignment.
        self.positions = self.run_numbering_algorithm_clustalo(
            base_sequence_id, alignment
        )
        logger.info(f"Positions computed: {self.positions}")

        # Create (or get) the StandardNumbering node in the database.
        StandardNumbering.get_or_save(
            name=f"{self.name}",
            definition=f"ClustalO based on base sequence {base_sequence_id}",
        )

        # Update the database with the relationships between proteins and standard numbering.
        self.save_positions(db)


if __name__ == "__main__":
    # Database connection setup parameters.
    uri = "bolt://129.69.129.130:7687"
    user = "neo4j"
    password = "12345678"

    from pyeed import Pyeed

    # Create a Pyeed object which automatically connects to the Neo4j database.
    eedb = Pyeed(uri, user, password)

    # Clear previous standard numbering relationships from the database.
    query = """
    MATCH (n:StandardNumbering)-[r:HAS_STANDARD_NUMBERING]-(c:Protein) DELETE r
    """
    eedb.db.execute_write(query)

    # Define sequences for testing the numbering algorithm.
    sequences = [
        ">seq1\nMTHKLLLTLLFTLLFSSAYSRG",
        ">seq2\nABCABCABCMTHKITLLLTLLFTLLFSSAYSRG",
        ">seq3\nMTHKILLLTLLFTLLFSSCYSRGARTHDB",
    ]

    proteins_dict = {
        "seq1": "MTHKLLLTLLFTLLFSSAYSRG",
        "seq2": "ABCABCABCMTHKITLLLTLLFTLLFSSAYSRG",
        "seq3": "MTHKILLLTLLFTLLFSSCYSRGARTHDB",
    }

    # Define a base sequence that will be used as the reference.
    base_sequence = {"id": "seq0", "sequence": "AMTHKLLLTLLFTLLFSSAYSRG"}

    from pyeed.tools.clustalo import ClustalOmega

    clustalO = ClustalOmega()

    # Insert the base sequence as the first sequence in the alignment list.
    sequences.insert(0, f">{base_sequence['id']}\n{base_sequence['sequence']}")

    # Create a dictionary for ClustalOmega from the sequences.
    sequences_dict = {
        base_sequence["id"]: base_sequence["sequence"],
        "seq1": "MTHKLLLTLLFTLLFSSAYSRG",
        "seq2": "ABCABCABCMTHKITLLLTLLFTLLFSSAYSRG",
        "seq3": "MTHKILLLTLLFTLLFSSCYSRGARTHDB",
    }

    # Perform multiple sequence alignment.
    alignment = clustalO.align(sequences_dict)

    # Instantiate the numbering tool and run the numbering algorithm.
    sn_tool = StandardNumberingTool("test_standard_numbering")
    sn_tool.positions = sn_tool.run_numbering_algorithm_pairwise("seq0", alignment)

    # Print a sample of the computed positions to verify the output.
    count = 0
    for i in sn_tool.positions:
        count += 1
        print(i, sn_tool.positions[i])
        # Only show the first few sequences
        if count > 10:
            break
