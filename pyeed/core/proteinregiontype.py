from enum import Enum


class ProteinRegionType(Enum):
    DOMAIN = "domain"
    SIGNAL_PEPTIDE = "signal peptide"
    TRANSMEMBRANE = "transmembrane"
    UNANNOTATED = "unannotated"
