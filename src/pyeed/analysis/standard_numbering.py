# here a standard numbering can be definied for proteins or dna sequences
# the defintion then points to all the dna or protein sequences that are part of the standard numbering
# on the relationship between the standard numbering and the dna or protein sequences, there is an array of positions

class StandardNumberingTool:
    def __init__(self, name):
        self.name = name
        self.positions = None

    def set_standard_numbering_with_given_base_sequence(self, base_sequence: dict[str, str], proteins_dict) -> dict[str, list[str]]:
        # base_sequence is a sequence which will be used to define the standard numbering
        # a clustalO alignment will be performed with the base_sequence and all the sequences in the proteins_dict
        # this will form as the basis for the standard numbering

        # the return will be a dictonary with the id as a key and the positions as a value

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


        return self.positions
    
        

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

    # the goal would be
    """
    
    """

    standard_numbering = StandardNumberingTool("test")
    positions = standard_numbering.set_standard_numbering_with_given_base_sequence(base_sequence, proteins_dict)

    for i in positions:
        print(i, positions[i])
