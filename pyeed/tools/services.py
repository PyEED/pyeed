from enum import Enum


class ServiceURL(Enum):
    CLUSTALO = "http://clustalo:5001/align"
    BLAST = "http://blast:6001/blast"
    FOLDSEEK = "http://foldseek:7001/foldseek"
    MMSEQS = "http://mmseqs:8001/easycluster"
