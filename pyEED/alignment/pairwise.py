from ctypes import alignment
from typing import List
from itertools import combinations
from Bio.Align import PairwiseAligner
from tqdm import tqdm

from pyEED.core import ProteinInfo
from pyEED.core import Alignment
from pyEED.core import StandardNumbering

from joblib import Parallel, delayed, cpu_count


def multi_pairwise_alignment(
    protein_infos: List[ProteinInfo],
    mode: str = "global",
    match: int = 1,
    mismatch: int = -1,
    gap_open: int = -1,
    gap_extend: int = 0,
    substitution_matrix: str = "None",
    n_jobs: int = None,
):
    pairs = list(combinations(protein_infos, 2))

    if n_jobs is None:
        n_jobs = cpu_count()

    alignments = Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(pairwise_alignment)(
            reference,
            query,
            mode,
            match,
            mismatch,
            gap_open,
            gap_extend,
            substitution_matrix,
        )
        for reference, query in tqdm(pairs, desc="⛓️ Aligning sequences")
    )

    return alignments


def pairwise_alignment(
    reference_info: ProteinInfo,
    query_info: ProteinInfo,
    mode: str = "global",
    match: int = 1,
    mismatch: int = -1,
    gap_open: int = -1,
    gap_extend: int = 0,
    substitution_matrix: str = "None",
) -> Alignment:
    """Creates a pairwise alignment between two sequences.

    Args:
        reference_info (ProteinInfo): Sequence 1 to be aligned.
        query_info (ProteinInfo): Sequence 2 to be aligned.
        mode (str, optional): If the alignment should be global or local. Defaults to "global".
        match (int, optional): Score of a match. Defaults to 1.
        mismatch (int, optional): Score of a mismatch. Defaults to -1.
        gap_open (int, optional): Gap open cost. Defaults to -1.
        gap_extend (int, optional): Gap extend cost. Defaults to 0.
        substitution_matrix (str, optional): Substitution matrix name. Defaults to "None".

    Raises:
        ValueError: If the mode is not global or local.
        FileNotFoundError: If the substitution matrix is not valid.

    Returns:
        Alignment: Pairwise alignment object with the alignment scores.
    """
    modes = ["global", "local"]

    if mode not in modes:
        raise ValueError(f"Invalid alignment mode: {mode}. Valid modes are: {modes}")

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend

    if substitution_matrix != "None":
        matrices = ["BLOSUM62", "BLOSUM45", "BLOSUM80", "PAM250", "PAM30", "PAM70"]
        try:
            from Bio.Align import substitution_matrices

            aligner.substitution_matrix = substitution_matrices.load(
                substitution_matrix
            )
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Invalid substitution matrix: {substitution_matrix}. Available matrices are: {matrices}"
            )

    alignment_result = aligner.align(reference_info.sequence, query_info.sequence)[0]
   
    #print(alignment_result)
    
    gaps = alignment_result.counts().gaps
    mismatches = alignment_result.counts().mismatches
    identities = alignment_result.counts().identities

    identity = identities / len(reference_info.sequence)

    standard_number = standard_numbering(
        reference_info=reference_info,
        alignment_result=alignment_result)
    
    alignment = Alignment(
        reference_seq=reference_info,
        query_seqs=[query_info],
        score=alignment_result.score,
        standard_numberings=[standard_number],
        identity=identity,
        gaps=gaps,
        mismatches=mismatches,
    )

    return alignment

def standard_numbering(
        reference_info: ProteinInfo,
        alignment_result: PairwiseAligner
        ) -> StandardNumbering:
    
    reference_seq_numbering = list(range(1, len(reference_info.sequence)+1))
    query_seq_numbering = []
    counter_gap = 1
    counter = 0

    for i in range(0,len(alignment_result[1][:])):
        if(alignment_result[0][i] == "-"):
            query_seq_numbering.append(str(reference_seq_numbering[counter]-1) + "." + str(counter_gap))
            counter_gap = counter_gap + 1
            
        elif(alignment_result[1][i] == "-"):
            counter_gap = 1
            counter = counter + 1
        
        else:
            query_seq_numbering.append(str(reference_seq_numbering[counter]))
            counter = counter + 1
            counter_gap = 1
  
    standard_number = StandardNumbering(
    sequence_id=reference_info.id,
    numbering=query_seq_numbering)
    
    return standard_number


def str_to_int(
        numbering: List
        ) -> List:
    
    reference_seq = []

    for number in numbering:

        find_point = number.find(".")

        if find_point != -1:
            reference_seq.append([int(number[0:find_point]), int(number[find_point+1:])])
        else:
            reference_seq.append([int(number),0])
    return reference_seq


def vizualize_standard_numbering(
        alignment: Alignment
    ):
    
    for j in range(len(alignment.query_seqs)):
    
        counter_ref = 0
        counter_query = 0

        identity = []
        print_target = []
        print_query = []

        reference_seq = list(range(1,len(alignment.reference_seq.sequence)+1))
        query_seq = str_to_int(alignment.standard_numberings[j].numbering)

        count_points = sum(1 for _, num in query_seq if num != 0)

        length = len(alignment.reference_seq.sequence) + count_points

        for i in range(length):

            if query_seq[counter_query][1] != 0:
                print_target.append('-')
                print_query.append(alignment.query_seqs[j].sequence[counter_query])
                counter_query += 1
            
            elif counter_ref != 0 and query_seq[counter_query][0] != reference_seq[counter_ref]:
                print_query.append('-')
                print_target.append(alignment.reference_seq.sequence[counter_ref])
                counter_ref += 1

            else: 
                print_query.append(alignment.query_seqs[j].sequence[counter_query])  
                print_target.append(alignment.reference_seq.sequence[counter_ref])
                counter_query += 1
                counter_ref += 1

        for i in range(0,length):  
            if print_target[i] == print_query[i]:
                identity.append('|')
            
            else:
                identity.append(' ')
    

    print(''.join(print_target))
    print(''.join(identity))
    print(''.join(print_query))                


                    
        

