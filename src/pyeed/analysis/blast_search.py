from typing import Optional

from Bio.Blast.Record import Blast
from Bio.Blast import NCBIWWW as BioBlastSearch

from pyeed.dbconnect import DatabaseConnector 
from pyeed.main import Pyeed 
from pyeed.tools.utility import chunks 
from rich.progess import Progress

#TODO: adapt description of the class
#TODO: MERGE in _to_db?? 
#TODO: set correct substitution matrix -> it is from the align module at the moment 

class BlastSearch:
    """
    Class for Blast search with biopython
    """
    
    def __init__(
        self,
        query: str,
        n_hits: int = 100,
        evalue: float = 1,
        substitution_matrix: str = "BLOSUM62",
        database: str = "nr",
        program: str = "blastp"         
    ) -> None:
        
        self.n_hits = n_hits
        self.evalue = evalue
        self.substitution_matrix = substitution_matrix
        self.database = database
        self.program = program
    
    def blast_search(
        self, 
        query_seq: str,
        db: Optional[DatabaseConnector] = None
    ) -> list[Blast]:
        """Performs a blast search with the provided query sequence.
        
        Args: 
            query_seq (str): The query sequence.
            db (Optional[DatabaseConnector], optional): A `DatabaseConnector` object. Defaults to None.
        
        Returns: 
            list[Blast]: List of blast results.
        """

        valid_programs = ["blastp", "blastn", "blastx", "tblastn", "tblastx"]
        valid_databases = ["nr", "nt", "swissprot", "pdb", "refseq_protein"]
        
        blast_search = self._get_blast_search()
        
        assert(
            self.program in valid_programs
        ), f"Invalid progam: {self.program}, valid programs are: {valid_programs}"
        
        assert(
            db in valid_databases
        ), f"Invalid database: {self.database}, valid databases are: {valid_databases}"
        
        result_handle = blast_search.qblast(sequence=query_seq)
        
        blast_records = list(NCBIXML.parse(result_handle))
        
        if db:
            self._to_db(blast_records, db)
        
        return blast_records
    
    def _to_db(
        self, 
        blast_records: list[Blast],
        db: DatabaseConnector
    ): 
        """Inserts the blast results into the database.

        Args:
            blast_records (list[Blast]): _description_
            db (DatabaseConnector): _description_
        """
        
        query = """
        UNWIND $blast_records AS record
        MATCH (p1:Protein {accession_id: record.query_id})
        MATCH (p2:Protein {accession_id: record.hit_id})
        
        """
        
        db.execute_write(query, {"blast_records": blast_records})
    
    def _get_blast_search(self) -> BioBlastSearch:
        
        blast_search = BioBlastSearch()
        blast_search.set_program(self.program)
        blast_search.set_matrix(self.substitution_matrix)
        blast_search.set_database(self.database)
        blast_search.set_expect(self.evalue)
        blast_search.set_hitlist_size(self.n_hits)
        
        if self.substitution_matrix != "None": 
            blast_search.substitution_matrix = self._load_substitution_matrix()
        
        return blast_search
    
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

    def _parse_blast_records(
        self, 
        blast_records: list[Blast],
    ) 

    def _load_substitution_matrix(self) -> "BioSubstitutionMatrix":
        from Bio.Align import substitution_matrices

        return substitution_matrices.load(self.substitution_matrix)
