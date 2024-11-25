from itertools import combinations
from typing import Dict, List, Optional

from Bio.Align import Alignment as Alignment
from Bio.Align import PairwiseAligner as BioPairwiseAligner
from Bio.Align.substitution_matrices import Array as BioSubstitutionMatrix
from joblib import Parallel, cpu_count, delayed
from pyeed.dbconnect import DatabaseConnector
from pyeed.main import Pyeed
from pyeed.tools.utility import chunks
from rich.progress import Progress


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

        results = aligner.align(list(seq1.values())[0], list(seq2.values())[0])

        return results[0]

    def align_pairwise(
        self,
        seq1: Dict[str, str],
        seq2: Dict[str, str],
        db: Optional[DatabaseConnector] = None,
    ) -> dict:
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
    ) -> Optional[List[dict]]:
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

        Returns:
            Optional[List[dict]]: A list of dictionaries containing the alignment results if
            `return_results` is True. If False, returns None.
        """

        # Fetch sequences if ids are provided
        if ids is not None and db is not None:
            sequences = self._get_id_sequence_dict(db, ids)

        if not sequences:
            raise ValueError(
                "Either sequences or ids (with a database connection) must be provided."
            )

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
                alignments = Parallel(n_jobs=cpu_count(), prefer="processes")(
                    delayed(self.align_pairwise)(
                        {pair[0]: sequences[pair[0]]},
                        {pair[1]: sequences[pair[1]]},
                        db=None,
                    )
                    for pair in pair_chunk
                )

                progress.update(align_task, advance=len(pair_chunk))

                if db:
                    self._to_db(alignments, db)
                    progress.update(db_task, advance=len(pair_chunk))

                if return_results:
                    all_alignments.extend(alignments)

        return all_alignments if return_results else None

    def _to_db(
        self,
        alignments: list[dict],
        db: DatabaseConnector,
    ):
        """Inserts the alignment results to pyeed graph database.

        Args:
            alignments (list[dict]): A list of dictionaries containing the alignment results.
            db (DatabaseConnector): A `DatabaseConnector` object.
        """

        query = """
        UNWIND $alignments AS alignment
        MATCH (p1:Protein {accession_id: alignment.query_id})
        MATCH (p2:Protein {accession_id: alignment.target_id})
        MERGE (p1)-[r:PAIRWISE_ALIGNED]->(p2)
        SET r.similarity = alignment.identity,
            r.mismatches = alignment.mismatches,
            r.gaps = alignment.gaps,
            r.score = alignment.score,
            r.query_aligned = alignment.query_aligned,
            r.target_aligned = alignment.target_aligned
        """

        db.execute_write(query, {"alignments": alignments})

    def _get_aligner(self) -> BioPairwiseAligner:
        """Creates a BioPython pairwise aligner object with the specified parameters
        from the class instance."""
        aligner = BioPairwiseAligner()
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
    ) -> dict:
        """Maps the alignment results to a dictionary.
        The dictionaly has the same signature as a `PairwiseAlignmentResult` object.

        Returns:
            dict: A dictionary containing the alignment results.
        """

        shorter_seq = min(alignment[0], alignment[1], key=lambda x: len(x))

        identities = alignment.counts().identities
        identity = identities / len(shorter_seq)
        gaps = alignment.counts().gaps
        mismatches = alignment.counts().mismatches

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

        if not ids:
            query = """
            MATCH (p:Protein)
            RETURN p.accession_id AS accession_id, p.sequence AS sequence
            """
            proteins = db.execute_read(query)
        else:
            query = """
            MATCH (p:Protein)
            WHERE p.accession_id IN $ids
            RETURN p.accession_id AS accession_id, p.sequence AS sequence
            """
            proteins = db.execute_read(query, {"ids": ids})

        return {protein["accession_id"]: protein["sequence"] for protein in proteins}

    def _load_substitution_matrix(self) -> "BioSubstitutionMatrix":
        from Bio.Align import substitution_matrices

        return substitution_matrices.load(self.substitution_matrix)


if __name__ == "__main__":
    s1 = "ACGTA"
    s2 = "ACGA"

    uri = "bolt://localhost:7688"
    username = "neo4j"
    password = "12345678"

    Pyeed(uri, user=username, password=password)
    
