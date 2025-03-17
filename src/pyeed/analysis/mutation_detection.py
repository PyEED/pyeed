from typing import Any, Optional

from loguru import logger
from pyeed.dbconnect import DatabaseConnector


class MutationDetection:
    def __init__(self) -> None:
        pass

    def get_sequence_data(
        self,
        sequence_id1: str,
        sequence_id2: str,
        db: DatabaseConnector,
        standard_numbering_tool_name: str,
        node_type: str = "Protein",
        region_based_sequence: Optional[str] = None,
    ) -> tuple[dict[str, str], dict[str, list[str]]]:
        """Fetch sequence and position data for two sequences from the database.

        Args:
            sequence_id1: First sequence accession ID
            sequence_id2: Second sequence accession ID
            db: Database connection instance
            standard_numbering_tool_name: Name of standard numbering tool to use
            node_type: Type of node to use (default: "Protein")
            region_based_sequence: Annotation of region to use for numbering (default: None)

        Returns:
            tuple containing:
                - dict[str, str]: Mapping of sequence IDs to sequences
                - dict[str, list[str]]: Mapping of sequence IDs to position lists

        Raises:
            ValueError: If standard numbering positions not found for both sequences
        """
        if region_based_sequence is not None:
            query = f"""
            MATCH (p:{node_type})-[rel:HAS_REGION]->(r:Region {{annotation: '{region_based_sequence}'}})-[rel2:HAS_STANDARD_NUMBERING]->(s:StandardNumbering)
            WHERE p.accession_id IN ['{sequence_id1}', '{sequence_id2}'] 
            AND s.name = '{standard_numbering_tool_name}'
            RETURN p.accession_id as id, p.sequence as sequence, rel2.positions as positions, rel.start as start, rel.end as end
            """
            results = db.execute_read(query)
        else:
            query = f"""
            MATCH (p:{node_type})-[r:HAS_STANDARD_NUMBERING]->(s:StandardNumbering)
            WHERE p.accession_id IN ['{sequence_id1}', '{sequence_id2}'] 
            AND s.name = '{standard_numbering_tool_name}'
            RETURN p.accession_id as id, p.sequence as sequence, r.positions as positions
            """
            results = db.execute_read(query)

        if len(results) < 2:
            raise ValueError(
                f"Could not find standard numbering positions for both sequences {sequence_id1} and {sequence_id2}"
            )
        if region_based_sequence is not None:
            sequences = {
                results[i]["id"]: results[i]["sequence"][
                    results[i]["start"] : results[i]["end"]
                ]
                for i in range(len(results))
            }
            positions = {
                results[i]["id"]: results[i]["positions"] for i in range(len(results))
            }

            return sequences, positions
        else:
            sequences = {result["id"]: result["sequence"] for result in results}
            positions = {result["id"]: result["positions"] for result in results}
            return sequences, positions

    def find_mutations(
        self,
        seq1: str,
        seq2: str,
        pos1: list[str],
        pos2: list[str],
    ) -> dict[str, Any]:
        """Compare two sequences and identify mutations between them.

        Args:
            seq1: First amino acid sequence
            seq2: Second amino acid sequence
            pos1: Standard numbering positions for first sequence
            pos2: Standard numbering positions for second sequence

        Returns:
            dict containing mutation information:
                - from_positions: List[int] - Source positions (1-based)
                - to_positions: List[int] - Target positions (1-based)
                - from_monomers: List[str] - Source amino acids
                - to_monomers: List[str] - Target amino acids
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

        return {
            "from_positions": from_positions,
            "to_positions": to_positions,
            "from_monomers": from_monomers,
            "to_monomers": to_monomers,
        }

    def save_mutations_to_db(
        self,
        mutations: dict[str, list[int | str]],
        db: DatabaseConnector,
        sequence_id1: str,
        sequence_id2: str,
        node_type: str = "Protein",
        region_based_sequence: Optional[str] = None,
    ) -> None:
        """Save detected mutations to the database.

        Args:
            mutations: Dictionary containing mutation information:
                - from_positions: List[int] - Source positions
                - to_positions: List[int] - Target positions
                - from_monomers: List[str] - Source amino acids
                - to_monomers: List[str] - Target amino acids
            db: Database connection instance
            sequence_id1: First sequence accession ID
            sequence_id2: Second sequence accession ID
            node_type: Type of node to use (default: "Protein")
            region_based_sequence: Annotation of region to use for numbering (default: None)
        """
        # Check if a mutation relationship already exists between these proteins
        if region_based_sequence is not None:
            query = f"""
            MATCH (p1:{node_type} {{accession_id: $sequence_id1}})-[rel:HAS_REGION]->(r1:Region {{annotation: $region_based_sequence}})-[rel_mutation:MUTATION]->(r2:Region {{annotation: $region_based_sequence}})<-[rel2:HAS_REGION]-(p2:{node_type} {{accession_id: $sequence_id2}})
            RETURN rel_mutation
            """
            existing_mutations = db.execute_read(
                query,
                {
                    "sequence_id1": sequence_id1,
                    "sequence_id2": sequence_id2,
                    "region_based_sequence": region_based_sequence,
                },
            )
        else:
            existing_mutations = db.execute_read(
                f"""
                MATCH (p1:{node_type})-[r:MUTATION]->(p2:{node_type})
                WHERE p1.accession_id = $sequence_id1 AND p2.accession_id = $sequence_id2
                RETURN r
                """,
                {"sequence_id1": sequence_id1, "sequence_id2": sequence_id2},
            )
        if existing_mutations:
            logger.debug(
                f"Mutation relationship already exists between {sequence_id1} and {sequence_id2}"
            )
            return

        if region_based_sequence is not None:
            # saving the mutation between the regions
            query = f"""
            MATCH (r1:Region {{annotation: $region_based_sequence}})<-[rel:HAS_REGION]-(p1:{node_type} {{accession_id: $sequence_id1}})
            MATCH (r2:Region {{annotation: $region_based_sequence}})<-[rel2:HAS_REGION]-(p2:{node_type} {{accession_id: $sequence_id2}})
            CREATE (r1)-[r:MUTATION]->(r2)
            SET r.from_positions = $from_positions,
                r.to_positions = $to_positions,
                r.from_monomers = $from_monomers,
                r.to_monomers = $to_monomers
            """
            params = {
                "sequence_id1": sequence_id1,
                "sequence_id2": sequence_id2,
                "region_based_sequence": region_based_sequence,
                "from_positions": mutations["from_positions"],
                "to_positions": mutations["to_positions"],
                "from_monomers": mutations["from_monomers"],
                "to_monomers": mutations["to_monomers"],
            }
            db.execute_write(query, params)
        else:
            query = f"""
            MATCH (p1:{node_type}), (p2:{node_type})
            WHERE p1.accession_id = $sequence_id1 AND p2.accession_id = $sequence_id2
            CREATE (p1)-[r:MUTATION]->(p2)
            SET r.from_positions = $from_positions,
                r.to_positions = $to_positions,
                r.from_monomers = $from_monomers,
                r.to_monomers = $to_monomers
            """
            params = {
                "sequence_id1": sequence_id1,
                "sequence_id2": sequence_id2,
                "from_positions": mutations["from_positions"],
                "to_positions": mutations["to_positions"],
                "from_monomers": mutations["from_monomers"],
                "to_monomers": mutations["to_monomers"],
            }
            db.execute_write(query, params)

        logger.debug(
            f"Saved {len(list(params['from_positions']))} mutations to database"
        )

    def get_mutations_between_sequences(
        self,
        sequence_id1: str,
        sequence_id2: str,
        db: DatabaseConnector,
        standard_numbering_tool_name: str,
        save_to_db: bool = True,
        debug: bool = False,
        node_type: str = "Protein",
        region_based_sequence: Optional[str] = None,
    ) -> dict[str, list[int | str]]:
        """Get mutations between two sequences using standard numbering.

        Args:
            sequence_id1: First sequence accession ID
            sequence_id2: Second sequence accession ID
            db: Database connection instance
            standard_numbering_tool_name: Name of standard numbering tool to use
            save_to_db: Whether to save mutations to database (default: True)
            node_type: Type of node to use (default: "Protein")
            region_based_sequence: Annotation of region to use for numbering (default: None)

        Returns:
            dict containing mutation information:
                - from_positions: List[int] - Source positions (1-based)
                - to_positions: List[int] - Target positions (1-based)
                - from_monomers: List[str] - Source amino acids
                - to_monomers: List[str] - Target amino acids

        Raises:
            ValueError: If standard numbering positions not found for both sequences
        """
        sequences, positions = self.get_sequence_data(
            sequence_id1,
            sequence_id2,
            db,
            standard_numbering_tool_name,
            node_type,
            region_based_sequence,
        )

        if debug:
            logger.info(f"Debug mode output: {sequences} and {positions}")

        mutations = self.find_mutations(
            sequences[sequence_id1],
            sequences[sequence_id2],
            positions[sequence_id1],
            positions[sequence_id2],
        )

        if save_to_db:
            self.save_mutations_to_db(
                mutations,
                db,
                sequence_id1,
                sequence_id2,
                node_type,
                region_based_sequence,
            )

        return mutations
