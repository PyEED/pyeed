# this tool should be able to detect the number of mutations between two sequences
# the two different sequences are either dna or protein and are already in the database
# how those two different sequences are found is not important for this tool, but should be done very carefully, since simply doing it between all of the possible sequences would be vastly out of scope

from pyeed.main import Pyeed
from pyeed.dbconnect import DatabaseConnector



class MutationDetection:
    def __init__(self):
        """
        Initialize the MutationDetection class.
        """
        None

    def get_mutations_between_sequences(self, sequence_id1: str, sequence_id2: str, db: DatabaseConnector, standard_numbering_tool_name: str, print_debug: bool = False, save_to_db: bool = True):
        """
        Get the mutations between two sequences using a standard numbering tool.

        Parameters:
            sequence_id1 (str): The accession ID of the first sequence
            sequence_id2 (str): The accession ID of the second sequence
            db (DatabaseConnector): The database connector object
            standard_numbering_tool_name (str): Name of the standard numbering tool to use
            print_debug (bool, optional): Whether to print debug information. Defaults to False.
            save_to_db (bool, optional): Whether to save mutations to database. Defaults to True.

        Returns:
            dict: A dictionary containing mutation information with the following structure:
                - sequence_id1: ID of the first sequence
                - sequence_id2: ID of the second sequence
                - from_positions: List of positions in the first sequence where mutations occur (1-based)
                - to_positions: List of positions in the second sequence where mutations occur (1-based)
                - from_monomers: List of original monomers at mutation positions
                - to_monomers: List of mutated monomers at mutation positions

        Raises:
            ValueError: If standard numbering positions cannot be found for both sequences

        Note:
            The mutations are detected by comparing sequence positions aligned using the 
            standard numbering tool. Only positions that exist in both sequences (common positions)
            are compared for mutations.
        """
        # Get the standard numbering positions for both sequences
        query = f"""
        MATCH (p:Protein)-[r:HAS_STANDARD_NUMBERING]->(s:StandardNumbering)
        WHERE p.accession_id IN ['{sequence_id1}', '{sequence_id2}'] 
        AND s.name = '{standard_numbering_tool_name}'
        RETURN p.accession_id as id, p.sequence as sequence, r.positions as positions
        """
        results = db.execute_read(query)

        if len(results) < 2:
            raise ValueError(f"Could not find standard numbering positions for both sequences {sequence_id1} and {sequence_id2} the results are: {results}")

        # Create dictionaries to store sequence and position info
        sequences = {}
        positions = {}
        for result in results:
            sequences[result['id']] = result['sequence']
            positions[result['id']] = result['positions']

        if print_debug:
            for i in sequences:
                print(f"{i}: {sequences[i]}")

            for i in positions:
                print(f"{i}: {positions[i]}")
        

        # Initialize lists to store mutation information
        from_positions = []
        to_positions = []
        from_monomers = []
        to_monomers = []

        # Compare sequences position by position using standard numbering
        seq1 = sequences[sequence_id1]
        seq2 = sequences[sequence_id2]
        pos1 = positions[sequence_id1]
        pos2 = positions[sequence_id2]

        # Create position-to-index mapping for both sequences
        pos_to_idx1 = {pos: idx for idx, pos in enumerate(pos1)}
        pos_to_idx2 = {pos: idx for idx, pos in enumerate(pos2)}

        # Find common positions between sequences
        common_positions = set(pos1) & set(pos2)

        # Check each common position for mutations
        for pos in common_positions:
            idx1 = pos_to_idx1[pos]
            idx2 = pos_to_idx2[pos]
            
            if seq1[idx1] != seq2[idx2]:
                from_positions.append(idx1 + 1)  # 1-based position
                to_positions.append(idx2 + 1)    # 1-based position
                from_monomers.append(seq1[idx1])
                to_monomers.append(seq2[idx2])

        data = {
            'sequence_id1': sequence_id1,
            'sequence_id2': sequence_id2,
            'from_positions': from_positions,
            'to_positions': to_positions,
            'from_monomers': from_monomers,
            'to_monomers': to_monomers
        }

        if save_to_db:
            self.save_mutations_to_db(data, db)

        return data
    

    def save_mutations_to_db(self, mutations: dict, db: DatabaseConnector):
        """
        Save detected mutations to the database as relationships between proteins.

        Parameters:
            mutations (dict): Dictionary containing mutation information with the structure:
                - sequence_id1: ID of the first sequence
                - sequence_id2: ID of the second sequence
                - from_positions: List of positions in the first sequence
                - to_positions: List of positions in the second sequence
                - from_monomers: List of original monomers
                - to_monomers: List of mutated monomers
            db (DatabaseConnector): The database connector object

        Note:
            Creates HAS_MUTATION relationships between proteins in the database,
            with properties storing the mutation details (positions and monomers).
        """

        # create the mutation relationship between the proteins
        for i in range(len(mutations['from_positions'])):
            query = f"""
            MATCH (p1:Protein), (p2:Protein)
            WHERE p1.accession_id = '{mutations['sequence_id1']}' AND p2.accession_id = '{mutations['sequence_id2']}'
            CREATE (p1)-[r:HAS_MUTATION]->(p2)
            SET r.from_position = {mutations['from_positions'][i]}, r.to_position = {mutations['to_positions'][i]}, r.from_monomer = '{mutations['from_monomers'][i]}', r.to_monomer = '{mutations['to_monomers'][i]}'
            """
            db.execute_write(query)
        


if __name__ == "__main__":
    uri = "bolt://localhost:7687"
    username = "neo4j"
    password = "12345678"

    file_path = "/home/nab/Niklas/TEM-lactamase/CARD_Data_Ontologies/aro.owl"

    eedb = Pyeed(uri, user=username, password=password)
    eedb.db.initialize_db_constraints(user=username, password=password)

    db = eedb.db

    name_of_standard_numbering_tool = 'test_standard_numbering'

    seq1 = 'KJO56189.1'
    seq2 = 'KLP91446.1'

    mutations = MutationDetection().get_mutations_between_sequences(seq1, seq2, db, name_of_standard_numbering_tool, print_debug=False)
    print(mutations)

    MutationDetection().save_mutations_to_db(mutations, db)
    
