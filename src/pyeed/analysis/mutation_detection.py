from loguru import logger
from pyeed.dbconnect import DatabaseConnector
from pyeed.model import Mutation


def get_sequence_data(
    sequence_id1: str,
    sequence_id2: str,
    db: DatabaseConnector,
    standard_numbering_tool_name: str,
) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Fetch sequence and position data for two sequences from the database.

    Args:
        sequence_id1: Accession ID of the first sequence
        sequence_id2: Accession ID of the second sequence
        db: Database connector instance
        standard_numbering_tool_name: Name of the standard numbering tool

    Returns:
        A tuple containing two dictionaries:
        - First dict mapping sequence IDs to sequences
        - Second dict mapping sequence IDs to position lists

    Raises:
        ValueError: If standard numbering positions are not found for both sequences
    """
    query = f"""
    MATCH (p:Protein)-[r:HAS_STANDARD_NUMBERING]->(s:StandardNumbering)
    WHERE p.accession_id IN ['{sequence_id1}', '{sequence_id2}'] 
    AND s.name = '{standard_numbering_tool_name}'
    RETURN p.accession_id as id, p.sequence as sequence, r.positions as positions
    """
    results = db.execute_read(query)

    if len(results) < 2:
        raise ValueError(
            f"Could not find standard numbering positions for both sequences {sequence_id1} and {sequence_id2}"
        )

    sequences = {result["id"]: result["sequence"] for result in results}
    positions = {result["id"]: result["positions"] for result in results}

    return sequences, positions


def find_mutations(
    seq1: str,
    seq2: str,
    pos1: list[str],
    pos2: list[str],
    sequence_id1: str,
    sequence_id2: str,
) -> Mutation:
    """Compare two sequences and find mutations between them.

    Args:
        seq1: First sequence string
        seq2: Second sequence string
        pos1: List of positions for first sequence
        pos2: List of positions for second sequence
        sequence_id1: Accession ID of the first sequence
        sequence_id2: Accession ID of the second sequence

    Returns:
        Mutation object containing the detected mutations
    """
    pos_to_idx1 = {pos: idx for idx, pos in enumerate(pos1)}
    pos_to_idx2 = {pos: idx for idx, pos in enumerate(pos2)}
    common_positions = set(pos1) & set(pos2)

    from_positions = []
    to_positions = []
    from_monomers = []
    to_monomers = []

    for pos in common_positions:
        idx1 = pos_to_idx1[pos]
        idx2 = pos_to_idx2[pos]

        if seq1[idx1] != seq2[idx2]:
            from_positions.append(idx1 + 1)  # 1-based position
            to_positions.append(idx2 + 1)  # 1-based position
            from_monomers.append(seq1[idx1])
            to_monomers.append(seq2[idx2])

    return Mutation(
        sequence_id1=sequence_id1,
        sequence_id2=sequence_id2,
        from_positions=from_positions,
        to_positions=to_positions,
        from_monomers=from_monomers,
        to_monomers=to_monomers,
    )


def save_mutations_to_db(mutations: Mutation, db: DatabaseConnector) -> None:
    """Save detected mutations to the database.

    Args:
        mutations: Mutation object containing mutation information
        db: Database connector instance
    """
    query = """
    MATCH (p1:Protein), (p2:Protein)
    WHERE p1.accession_id = $sequence_id1 AND p2.accession_id = $sequence_id2
    CREATE (p1)-[r:MUTATION]->(p2)
    SET r.from_positions = $from_positions,
        r.to_positions = $to_positions,
        r.from_monomers = $from_monomers,
        r.to_monomers = $to_monomers
    """
    params = {
        "sequence_id1": mutations.sequence_id1,
        "sequence_id2": mutations.sequence_id2,
        "from_positions": mutations.from_positions,
        "to_positions": mutations.to_positions,
        "from_monomers": mutations.from_monomers,
        "to_monomers": mutations.to_monomers,
    }
    db.execute_write(query, params)
    logger.debug(f"Saved {len(mutations.from_positions)} mutations to database")


def get_mutations_between_sequences(
    sequence_id1: str,
    sequence_id2: str,
    db: DatabaseConnector,
    standard_numbering_tool_name: str,
    save_to_db: bool = True,
) -> Mutation:
    """Get the mutations between two sequences using a standard numbering tool.

    Args:
        sequence_id1: Accession ID of the first sequence
        sequence_id2: Accession ID of the second sequence
        db: Database connector instance
        standard_numbering_tool_name: Name of the standard numbering tool to use
        save_to_db: Whether to save mutations to database

    Returns:
        Mutation object containing the detected mutations
    """
    sequences, positions = get_sequence_data(
        sequence_id1, sequence_id2, db, standard_numbering_tool_name
    )

    mutations = find_mutations(
        sequences[sequence_id1],
        sequences[sequence_id2],
        positions[sequence_id1],
        positions[sequence_id2],
        sequence_id1,
        sequence_id2,
    )

    if save_to_db:
        save_mutations_to_db(mutations, db)

    return mutations


if __name__ == "__main__":
    from pyeed.main import Pyeed

    uri = "bolt://localhost:7687"
    username = "neo4j"
    password = "12345678"

    eedb = Pyeed(uri, user=username, password=password)
    eedb.db.initialize_db_constraints(user=username, password=password)

    seq1 = "KJO56189.1"
    seq2 = "KLP91446.1"
    name_of_standard_numbering_tool = "test_standard_numbering"

    mutations = get_mutations_between_sequences(
        seq1, seq2, eedb.db, name_of_standard_numbering_tool
    )
    logger.info(f"Found mutations: {mutations}")
