# here a standard numbering can be definied for proteins or dna sequences
# the defintion then points to all the dna or protein sequences that are part of the standard numbering
# on the relationship between the standard numbering and the dna or protein sequences, there is an array of positions
from loguru import logger

from pyeed.model import StandardNumbering
from pyeed.dbconnect import DatabaseConnector


class StandardNumberingTool:
    def __init__(self, name):
        self.name = name
        self.positions = None

    def _get_proteins_dict(self, db: DatabaseConnector):
        """
        This function will get all the proteins from the database and return them as a dictionary
        """
        query = """
        MATCH (p:Protein)
        WHERE p.sequence IS NOT NULL
        RETURN p.accession_id AS accession_id, p.sequence AS sequence
        """
        proteins_read = db.execute_read(query)
        proteins_dict = {
            protein["accession_id"]: protein["sequence"] for protein in proteins_read
        }

        return proteins_dict

    def get_protein_base_sequence(self, base_sequence_id: str, db: DatabaseConnector):
        """
        This function will get the base sequence from the database
        """
        query = f"""
        MATCH (p:Protein)
        WHERE p.accession_id = '{base_sequence_id}'
        RETURN p.accession_id AS accession_id, p.sequence AS sequence
        """
        base_sequence_read = db.execute_read(query)
        base_sequence = {
            "id": base_sequence_read[0]["accession_id"],
            "sequence": base_sequence_read[0]["sequence"]
        }

        return base_sequence
    
    def safe_positions(self, db: DatabaseConnector):
        """
        This function will save the positions to the database
        The positions need to be one long string
        """
        for protein_id in self.positions:
            query = f"""
                MATCH (p:Protein {{accession_id: '{protein_id}'}})
                MATCH (s:StandardNumbering {{name: '{self.name}'}})
                MERGE (p)-[r:HAS_STANDARD_NUMBERING]->(s)
                SET r.positions = {str(self.positions[protein_id])}
            """
            db.execute_write(query)

    def run_numbering_algorithm(self, base_sequence_id: str, alignment):
        """
        Core algorithm to compute numbering positions for all sequences relative to the base sequence.
        """
        logger.info(f"Running numbering algorithm for base sequence {base_sequence_id}")

        positions = {}

        # the base sequence has the postions 1, ...
        base_sequence = alignment[0].seq
        positions[base_sequence_id] = [str(i + 1) for i in range(len(base_sequence.replace('-', '')))]

        base_seq_counter = -1
        # iterate over the positions in the alignment
        for pos in range(len(alignment[0].seq)):

            # this is the amino acid in the base sequence
            base_aa = alignment[0].seq[pos]

            if base_aa != '-':
                base_seq_counter += 1
            
            for j in range(1, len(alignment)):
                # get the sequence
                sequence = alignment[j].seq
                # get the id
                sequence_id = alignment[j].id

                if base_aa == '-' and sequence[pos] == '-':
                    # both are inserts
                    continue

                if base_aa == '-' and sequence[pos] != '-':
                    # insert in the base sequence
                    if sequence_id not in positions:
                        positions[sequence_id] = []
                        positions[sequence_id].append('0.0')
                    
                    else:
                        # check if the previous position is an insert
                        if '.' in positions[sequence_id][-1]:
                            # get the number of the insert
                            insert_number = int(positions[sequence_id][-1].split('.')[1])
                            base_number = int(positions[sequence_id][-1].split('.')[0])
                            positions[sequence_id].append(f"{base_number}.{insert_number + 1}")
                        
                        else:
                            positions[sequence_id].append(f"{positions[sequence_id][-1]}.1")

                if base_aa != '-' and sequence[pos] == '-':
                    if sequence_id not in positions:
                        positions[sequence_id] = [positions[base_sequence_id][base_seq_counter] + '.GAP']
                        
                    positions[sequence_id].append(f"{positions[sequence_id][-1].split('.')[0]}.GAP")

                if base_aa != '-' and sequence[pos] != '-':
                    # both are amino acids
                    if sequence_id not in positions:
                        positions[sequence_id] = []
                        positions[sequence_id].append(positions[base_sequence_id][base_seq_counter])
                    
                    else:
                        # check the previous position base number (does not matter if it is an insert)
                        # increment the number
                        if '.' in positions[sequence_id][-1] and 'GAP' not in positions[sequence_id][-1]:
                            # get the number of the insert
                            base_number = int(positions[sequence_id][-1].split('.')[0])
                            positions[sequence_id].append(f"{base_number + 1}")
                        
                        else:
                            positions[sequence_id].append(f"{int(positions[base_sequence_id][base_seq_counter])}")


        # clean up all the GAP and remove them
        for protein_id in positions:
            positions[protein_id] = [pos for pos in positions[protein_id] if 'GAP' not in pos]

                
        return positions




    def apply_standard_numbering(self, base_sequence_id: str, db: DatabaseConnector):
        """
        This function will set the standard numbering for a given base sequence and all the sequences in the database
        The sequences will be aligned with clustal omega and the positions will be determined

        base_sequence_id: str -> the id of the base sequence
        db: DatabaseConnector -> the database connector object

        """
        # get all the proteins from the database
        proteins_dict = self._get_proteins_dict(db)

        # get the base sequence
        base_sequence = self.get_protein_base_sequence(base_sequence_id, db)
        print(f"Base sequence: {base_sequence}")

        # remove the base sequence from the proteins_dict
        proteins_dict.pop(base_sequence_id)

        
        print(f"Base sequence: {base_sequence}")
        
        # run clustal omega
        from pyeed.tools.clustalo import ClustalOmega

        # create a list of sequences
        sequences = []
        sequences.append(f">{base_sequence['id']}\n{base_sequence['sequence']}")
        for key in proteins_dict:
            # format is >id\nsequence
            sequences.append(f">{key}\n{proteins_dict[key]}")

        print(sequences[:3])
        
        # run clustal omega
        clustalO = ClustalOmega()
        alignment = clustalO.align(sequences)

        # get the alignment of base sequence is the first sequence
        base_sequence_alignment = alignment[0].seq

        # get all postions relative to the base sequence
        positions_base = list(range(len(base_sequence_alignment.replace('-', ''))))
        positions_base = [str(i) for i in positions_base]
        
        # get all positions for all the sequences
        # we want to use inserts as 2.1 .. to check wether it is an insert we take a look at the base sequence at the postions if - is present
        self.positions = self.run_numbering_algorithm(base_sequence_id, alignment)

        print(self.positions['AEQ59620.1'])
        print(list(self.positions.values())[:2])


        # update the database with the standard numbering
        # create the standard numbering node
        StandardNumbering.get_or_save(
            name=f"{self.name}", definition=f"ClustalO based on base sequence {base_sequence_id}"
        )

        # create the relationships between the standard numbering and the proteins
        self.safe_positions(db)
    
        

if __name__ == "__main__":


    uri = "bolt://127.0.0.1:7687"
    user = "neo4j"
    password = "12345678"

    from pyeed import Pyeed


    # Create a Pyeed object, automatically connecting to the database
    eedb = Pyeed(uri, user, password)

    # remove all relationships on the previous standard numbering
    query = """
    MATCH (n:StandardNumbering)-[r:HAS_STANDARD_NUMBERING]-(c:Protein) DELETE r
    """

    eedb.db.execute_write(query)

    
    # test the standard numbering
    sequences = [
        ">seq1\nMTHKLLLTLLFTLLFSSAYSRG",
        ">seq2\nABCABCABCMTHKITLLLTLLFTLLFSSAYSRG",
        ">seq3\nMTHKILLLTLLFTLLFSSCYSRGARTHDB",
    ]

    proteins_dict = {
        "seq1": "MTHKLLLTLLFTLLFSSAYSRG",
        "seq2": "ABCABCABCMTHKITLLLTLLFTLLFSSAYSRG",
        "seq3": "MTHKILLLTLLFTLLFSSCYSRGARTHDB"
    }

    base_sequence = {
        "id": 'seq0',
        "sequence": 'AMTHKLLLTLLFTLLFSSAYSRG'
    }

    from pyeed.tools.clustalo import ClustalOmega

    clustalO = ClustalOmega()

    # add base sequence to the sequences
    sequences.insert(0, f">{base_sequence['id']}\n{base_sequence['sequence']}")

    alignment = clustalO.align(sequences)

    sn_tool = StandardNumberingTool("test_standard_numbering")
    sn_tool.positions = sn_tool.run_numbering_algorithm('seq0', alignment)    

    count = 0
    for i in sn_tool.positions:
        count += 1
        print(i, sn_tool.positions[i])

        if count > 10:
            break
    

