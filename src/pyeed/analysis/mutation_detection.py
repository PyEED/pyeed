# this tool should be able to detect the number of mutations between two sequences
# the two different sequences are either dna or protein and are already in the database
# how those two different sequences are found is not important for this tool, but should be done very carefully, since simply doing it between all of the possible sequences would be vastly out of scope

from pyeed.main import Pyeed
from pyeed.dbconnect import DatabaseConnector



class MutationDetection:
    def __init__(self):
        None

    def get_number_of_mutations(self, sequence_id1: str, sequence_id2: str, db: DatabaseConnector):
        """
        This function will get the number of mutations between two sequences.
        
        Parameters:
        sequence_id1 (str): The sequence id for the first sequence
        sequence_id2 (str): The sequence id for the second sequence
        db (DatabaseConnector): The database connector object

        Returns:
        int: The number of mutations
        """

        # get the sequence for the first sequence id
        query = f"""
        MATCH (p:Protein {{accession_id: '{sequence_id1}'}})
        RETURN p.sequence AS sequence
        """
        sequence_read = db.execute_read(query)
        sequence1 = sequence_read[0]['sequence']

        # get the sequence for the second sequence id
        query = f"""
        MATCH (p:Protein {{accession_id: '{sequence_id2}'}})
        RETURN p.sequence AS sequence
        """
        sequence_read = db.execute_read(query)
        sequence2 = sequence_read[0]['sequence']

        # get the number of mutations
        number_of_mutations = 0
        for i in range(len(sequence1)):
            if sequence1[i] != sequence2[i]:
                number_of_mutations += 1

        return number_of_mutations
    

if __name__ == "__main__":
    uri = "bolt://localhost:7687"
    username = "neo4j"
    password = "12345678"

    file_path = "/home/nab/Niklas/TEM-lactamase/CARD_Data_Ontologies/aro.owl"

    eedb = Pyeed(uri, user=username, password=password)
    eedb.db.initialize_db_constraints(user=username, password=password)

    db = eedb.db
    