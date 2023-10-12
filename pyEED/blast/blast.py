from typing import List
from Bio.Blast import NCBIWWW, NCBIXML

from pyeed.core.proteinsequence import ProteinSequence


def pblast(sequence: str, n_hits: int = None):
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence, hitlist_size=n_hits)
    blast_record = NCBIXML.read(result_handle)

    accessions = _get_accessions(blast_record)

    sequences = []
    for acc in accessions:
        sequences.append(ProteinSequence.from_ncbi(acc))

    return sequences


def nblast(sequence: str, n_hits: int = None) -> List[ProteinSequence]:
    result_handle = NCBIWWW.qblast("blastn", "nr", sequence, hitlist_size=n_hits)
    blast_record = NCBIXML.read(result_handle)


def _get_accessions(blast_record: NCBIXML) -> List[str]:
    accessions = []
    for alignment in blast_record.alignments:
        accessions.append(alignment.accession)
    return accessions
