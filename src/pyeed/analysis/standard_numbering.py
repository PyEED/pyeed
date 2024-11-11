# here a standard numbering can be definied for proteins or dna sequences
# the defintion then points to all the dna or protein sequences that are part of the standard numbering
# on the relationship between the standard numbering and the dna or protein sequences, there is an array of positions
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

    def _get_protein_base_sequence(self, base_sequence_id: str, db: DatabaseConnector):
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
    
    def _safe_positions(self, db: DatabaseConnector):
        """
        This function will save the positions to the database
        """
        for protein_id in self.positions:
            query = f"""
                MATCH (p:Protein {{accession_id: '{protein_id}'}})
                MATCH (s:StandardNumbering {{name: '{self.name}'}})
                MERGE (p)-[r:HAS_STANDARD_NUMBERING]->(s)
                SET r.positions = {self.positions[protein_id]}
            """
            db.execute_write(query)

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
        base_sequence = self._get_protein_base_sequence(base_sequence_id, db)
        
        # run clustal omega
        from pyeed.tools.clustalo import ClustalOmega

        # create a list of sequences
        sequences = []
        sequences.append(f">{base_sequence['id']}\n{base_sequence['sequence']}")
        for key in proteins_dict:
            # format is >id\nsequence
            sequences.append(f">{key}\n{proteins_dict[key]}")
        
        # run clustal omega
        clustalO = ClustalOmega()
        alignment = clustalO.align(sequences)

        # get the alignment of base sequence is the first sequence
        base_sequence_alignment = alignment[0].seq

        # get all postions relative to the base sequence
        positions_base = list(range(len(base_sequence_alignment.replace('-', ''))))
        positions_base = [i for i in positions_base]
        
        # get all positions for all the sequences
        # we want to use inserts as 2.1 .. to check wether it is an insert we take a look at the base sequence at the postions if - is present
        positions = {}
        for i in range(1, len(alignment)):
            sequence = alignment[i].seq
            positions_sequence: list[str] = []

            for j in range(len(base_sequence_alignment.replace('-', ''))):
                if base_sequence_alignment[j] == '-' and sequence[j] != '-':
                    # insert, but could be first second etc insert
                    if len(positions_sequence) > 0:
                        if positions_sequence[-1].count('.') != 0:
                            if positions_sequence[-1].split('.')[0] == str(j):
                                insert_number = int(positions_sequence[-1].split('.')[1]) + 1
                            else:
                                insert_number = 1
                        else:
                            insert_number = 1
                    else:
                        insert_number = 1
                        
                    for k in range(j, len(sequence)):
                        if base_sequence_alignment[k] != '-':
                            break
                        else:
                            positions_sequence.append(str(j) + '.' + str(insert_number+(k-j)))
                else:
                    positions_sequence.append(str(j))


            positions[alignment[i].id] = positions_sequence
        
        

        positions[base_sequence['id']] = [str(i) for i in positions_base]

        self.positions = positions


        # update the database with the standard numbering
        # create the standard numbering node
        StandardNumbering.get_or_save(
            name=f"{self.name}", definition=f"ClustalO based on base sequence {base_sequence_id}"
        )

        # create the relationships between the standard numbering and the proteins
        self._safe_positions(db)
    
        

if __name__ == "__main__":

    # test the standard numbering
    sequences = [
        ">seq1\nMTHKLLLTLLFTLLFSSAYSRG",
        ">seq2\nMTHKITLLLTLLFTLLFSSAYSRG",
        ">seq3\nMTHKILLLTLLFTLLFSSCYSRG",
    ]

    base_sequence = {
        "id": "seq0",
        "sequence": "MTHKLLLTLLFTLLFSSTYSRG"
    }

    proteins_dict = {
        "seq1": "MTHKLLLTLLFTLLFSSAYSRG",
        "seq2": "MTHKITLLLTLLFTLLFSSAYSRG",
        "seq3": "MTHKILLLTLLFTLLFSSCYSRG"
    }

    None
