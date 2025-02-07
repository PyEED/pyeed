from enum import Enum


class ServiceURL(Enum):
    CLUSTALO = "http://129.69.129.130:5001/align"
    BLAST = "http://blast:6001/"
    FOLDSEEK = "http://foldseek:7001/foldseek"
    MMSEQS = "http://mmseqs:8001/easycluster"
