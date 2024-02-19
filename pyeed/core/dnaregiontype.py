from enum import Enum


class DNARegionType(Enum):
    CODING_SEQUENCE = "coding sequence"
    EXON = "exon"
    INTRON = "intron"
    GENE = "gene"
    PROMOTER = "promoter"
    ENHANCER = "enhancer"
    UNANNOTATED = "unannotated"
