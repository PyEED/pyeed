from itertools import combinations
from typing import Any, Dict, Optional

from Bio.Align import Alignment as Alignment
from Bio.Align import PairwiseAligner as BioPairwiseAligner
from Bio.Align.substitution_matrices import Array as BioSubstitutionMatrix
from joblib import Parallel, cpu_count, delayed
from loguru import logger
from rich.progress import Progress

from pyeed.dbconnect import DatabaseConnector
from pyeed.tools.utility import chunks


class PairwiseAligner:
    """A class for pairwise and parallelized multi-pairwise sequence
    alignments using BioPython. Align methods can write the results directly
    to a neo4j database.
    """

    def __init__(
        self,
        mode: str = "global",
        match: int = 1,
        mismatch: int = -1,
        gap_open: int = -1,
        gap_exted: int = 0,
        substitution_matrix: str = "None",
    ) -> None:
        self.mode = mode
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extend = gap_exted
        self.substitution_matrix = substitution_matrix

    def _align(
        self,
        seq1: Dict[str, str],
        seq2: Dict[str, str],
    ) -> Alignment:
        """Aligns two sequences and returns the alignment results.

        Args:
            seq1 (Dict[str, str]): Sequence 1 to align. Key is the sequence ID.
            seq2 (Dict[str, str]): Sequence 2 to align. Key is the sequence ID.
            progress (Progress, optional): `Rich` progress. Defaults to None.
            task_id (TaskID, optional): `Rich` task_id. Defaults to None.

        Returns:
            Alignment: The alignment results.
        """

        aligner = self._get_aligner()

        results = aligner.align(list(seq1.values())[0], list(seq2.values())[0])  # type: ignore

        return results[0]  # type: ignore

    def align_pairwise(
        self,
        seq1: Dict[str, str],
        seq2: Dict[str, str],
        db: Optional[DatabaseConnector] = None,
    ) -> dict[str, Any]:
        """Aligns two sequences and returns the alignment results.
        If a `DatabaseConnector` object is provided, the results are added to the database.

        Args:
            seq1 (Dict[str, str]): Sequence 1 to align. Key is the sequence ID.
            seq2 (Dict[str, str]): Sequence 2 to align. Key is the sequence ID.
            db (DatabaseConnector, optional): A `DatabaseConnector` object. Defaults to None.

        Returns:
            dict: Has the same signature as a `PairwiseAlignmentResult` object.
        """

        alignment_result = self._align(seq1, seq2)

        alignment = self._map_alignment_results(alignment_result, seq1, seq2)

        if db:
            self._to_db([alignment], db)

        return alignment

    def align_multipairwise(
        self,
        ids: Optional[list[str]] = None,
        sequences: Optional[dict[str, str]] = None,
        db: Optional[DatabaseConnector] = None,
        batch_size: int = 500,
        return_results: bool = True,
        pairs: Optional[list[tuple[str, str]]] = None,
        node_type: str = "Protein",
        region_ids_neo4j: Optional[list[str]] = None,
        num_cores: int = cpu_count() - 1,
    ) -> Optional[list[dict[str, Any]]]:
        """
        Creates all possible pairwise alignments from a dictionary of sequences or from sequence IDs.
        If a `DatabaseConnector` object is provided, the results can be added to the database.

        If `ids` are provided, the method fetches sequences from the database using these IDs.
        If neither `sequences` nor `ids` are provided, a `ValueError` is raised.

        Args:
            ids (Optional[list[str]]): A list of sequence IDs to fetch and align. If provided,
                sequences will be retrieved from the database if a `DatabaseConnector` is provided.
            sequences (Optional[dict[str, str]]): A dictionary of sequences to align, where the key
                is the sequence ID and the value is the sequence string. If provided, these sequences
                are used directly for alignment.
            db (Optional[DatabaseConnector]): A `DatabaseConnector` object for fetching sequences
                if `ids` are provided. Defaults to None.
            batch_size (int): The size of the chunks to process for alignments. Defaults to 500.
            return_results (bool): Whether to return the list of alignment results. If set to False,
                the function will perform the alignments and any database insertion without
                returning the results, which can reduce memory usage. Defaults to True.
            pairs (Optional[list[tuple[str, str]]]): A list of tuples, where each tuple contains two
                sequence IDs to align. If provided, only these pairs will be aligned.
            node_type (str): The type of node to align. Defaults to "Protein".
            region_ids_neo4j (Optional[list[str]]): A list of region IDs for the sequence cuting based on region_based_sequence.
        Returns:
            Optional[List[dict]]: A list of dictionaries containing the alignment results if
            `return_results` is True. If False, returns None.
        """

        # Fetch sequences if ids are provided
        if ids is not None and db is not None:
            sequences = self._get_id_sequence_dict(db, ids, node_type, region_ids_neo4j)

        logger.info(
            f"Length of sequences: {len(sequences)} and length of pairs: {len(pairs)} and length of ids: {len(ids)} and the length of the region_ids_neo4j: {len(region_ids_neo4j)}"
        )
        logger.info(f"IDS: {ids}")
        logger.info(f"Region IDs: {region_ids_neo4j}")
        logger.info(f"Pairs: {pairs}")
        logger.info(f"Sequences: {sequences.keys()}")

        if not sequences:
            raise ValueError(
                "Either sequences or ids (with a database connection) must be provided."
            )

        if pairs is None:
            pairs = list(combinations(sequences, 2))

        total_pairs = len(pairs)
        all_alignments = []

        with Progress() as progress:
            align_task = progress.add_task(
                f"⛓️ Aligning {total_pairs} sequence pairs...", total=total_pairs
            )
            db_task = progress.add_task(
                "📥 Inserting alignment results to database...", total=total_pairs
            )

            for pair_chunk in chunks(pairs, batch_size):
                # Align the pairs in the current chunk
                alignments = Parallel(n_jobs=num_cores, prefer="processes")(
                    delayed(self.align_pairwise)(
                        {pair[0]: sequences[pair[0]]},
                        {pair[1]: sequences[pair[1]]},
                        db=None,
                    )
                    for pair in pair_chunk
                )

                progress.update(align_task, advance=len(pair_chunk))

                if db:
                    self._to_db(alignments, db, node_type, region_ids_neo4j)
                    progress.update(db_task, advance=len(pair_chunk))

                if return_results:
                    all_alignments.extend(alignments)

        return all_alignments if return_results else None

    def _to_db(
        self,
        alignments: list[dict[str, Any]],
        db: DatabaseConnector,
        node_type: str = "Protein",
        region_ids_neo4j: Optional[list[str]] = None,
    ) -> None:
        """Inserts the alignment results to pyeed graph database.

        Args:
            alignments (list[dict]): A list of dictionaries containing the alignment results.
            db (DatabaseConnector): A `DatabaseConnector` object.
            node_type (str): The type of node to align. Defaults to "Protein".
            region_ids_neo4j (Optional[list[str]]): A list of region IDs for the sequence cuting based on region_based_sequence.
        """

        if region_ids_neo4j is None:
            query = f"""
            UNWIND $alignments AS alignment
            MATCH (p1:{node_type} {{accession_id: alignment.query_id}})
            MATCH (p2:{node_type} {{accession_id: alignment.target_id}})
            MERGE (p1)-[r:PAIRWISE_ALIGNED]->(p2)
            SET r.similarity = alignment.identity,
            r.mismatches = alignment.mismatches,
            r.gaps = alignment.gaps,
            r.score = alignment.score,
            r.query_aligned = alignment.query_aligned,
            r.target_aligned = alignment.target_aligned
            """
            db.execute_write(query, parameters={"alignments": alignments})
        else:
            query = f"""
            UNWIND $alignments AS alignment
            MATCH (p1:{node_type} {{accession_id: alignment.query_id}})-[rel1:HAS_REGION]->(r1:Region)
            MATCH (p2:{node_type} {{accession_id: alignment.target_id}})-[rel2:HAS_REGION]->(r2:Region)
            WHERE id(r1) IN $region_ids_neo4j AND id(r2) IN $region_ids_neo4j
            MERGE (r1)-[r:PAIRWISE_ALIGNED]->(r2)
            SET r.similarity = alignment.identity,
            r.mismatches = alignment.mismatches,
            r.gaps = alignment.gaps,
            r.score = alignment.score,
            r.query_aligned = alignment.query_aligned,
            r.target_aligned = alignment.target_aligned
            """
            db.execute_write(
                query,
                parameters={
                    "alignments": alignments,
                    "region_ids_neo4j": region_ids_neo4j,
                },
            )

    def _get_aligner(self) -> BioPairwiseAligner:
        """Creates a BioPython pairwise aligner object with the specified parameters
        from the class instance."""
        aligner = BioPairwiseAligner()  # type: ignore
        aligner.mode = self.mode
        aligner.match_score = self.match
        aligner.mismatch_score = self.mismatch
        aligner.open_gap_score = self.gap_open
        aligner.extend_gap_score = self.gap_extend

        if self.substitution_matrix != "None":
            aligner.substitution_matrix = self._load_substitution_matrix()

        return aligner

    def _map_alignment_results(
        self, alignment: Alignment, seq1: Dict[str, str], seq2: Dict[str, str]
    ) -> dict[str, Any]:
        """Maps the alignment results to a dictionary.
        The dictionaly has the same signature as a `PairwiseAlignmentResult` object.

        Returns:
            dict: A dictionary containing the alignment results.
        """

        shorter_seq = min(alignment[0], alignment[1], key=lambda x: len(x))

        identities = alignment.counts().identities  # type: ignore
        identity = identities / len(shorter_seq)
        gaps = alignment.counts().gaps  # type: ignore
        mismatches = alignment.counts().mismatches  # type: ignore

        query_aligned = alignment[0]
        target_aligned = alignment[1]

        result_dict = {
            "query_id": list(seq1.keys())[0],
            "target_id": list(seq2.keys())[0],
            "score": alignment.score,  # type: ignore
            "identity": identity,
            "gaps": gaps,
            "mismatches": mismatches,
            "query_aligned": query_aligned,
            "target_aligned": target_aligned,
        }

        return result_dict

    def _get_id_sequence_dict(
        self,
        db: DatabaseConnector,
        ids: list[str] = [],
        node_type: str = "Protein",
        region_ids_neo4j: Optional[list[str]] = None,
    ) -> dict[str, str]:
        """Gets all sequences from the database and returns them in a dictionary.
        Key is the accession id and value is the sequence.
        If no ids are provided, all sequences are returned.

        Args:
            db (DatabaseConnector): A `DatabaseConnector` object.
            ids (list[str], optional): List of accession ids to fetch. Defaults to [].

        Returns:
            dict[str, str]: Dictionary of sequences with accession id as key.
        """

        if ids != []:
            if region_ids_neo4j is not None:
                query = f"""
                MATCH (p:{node_type})-[e:HAS_REGION]->(r:Region)
                WHERE id(r) IN {region_ids_neo4j} AND p.accession_id IN {ids}
                RETURN p.accession_id AS accession_id, e.start AS start, e.end AS end, p.sequence AS sequence
                """
                nodes = db.execute_read(query)

            else:
                query = f"""
                MATCH (p:{node_type})
                WHERE p.accession_id IN $ids
                RETURN p.accession_id AS accession_id, p.sequence AS sequence
                """
                nodes = db.execute_read(query, parameters={"ids": ids})

        else:
            if region_ids_neo4j is not None:
                query = f"""
                MATCH (p:{node_type})-[e:HAS_REGION]->(r:Region)
                WHERE id(r) IN $region_ids_neo4j
                RETURN p.accession_id AS accession_id, e.start AS start, e.end AS end, p.sequence AS sequence
                """
                nodes = db.execute_read(
                    query,
                    parameters={
                        "region_ids_neo4j": region_ids_neo4j,
                    },
                )
            else:
                query = f"""
                MATCH (p:{node_type})
                RETURN p.accession_id AS accession_id, p.sequence AS sequence
                """
                nodes = db.execute_read(query)

        if region_ids_neo4j is not None:
            return {
                node["accession_id"]: node["sequence"][node["start"] : node["end"]]
                for node in nodes
            }
        else:
            return {node["accession_id"]: node["sequence"] for node in nodes}

    def _load_substitution_matrix(self) -> "BioSubstitutionMatrix":
        from Bio.Align import substitution_matrices

        return substitution_matrices.load(self.substitution_matrix)  # type: ignore
