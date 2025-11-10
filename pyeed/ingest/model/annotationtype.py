from enum import Enum


class AnnotationType(str, Enum):
    """Protein and DNA annotation types."""

    SITE = "site"
    ACTIVE_SITE = "active_site"
    ALLOSTERIC_SITE = "allosteric_site"
    ALPHAHELIX = "alpha_helix"
    BETASTRAND = "beta_strand"
    BINDING_SITE = "binding_site"
    MATURE_PROTEIN = "mature_protein"
    CODING_SEQ = "coding_sequence"
    DNA = "DNA"
    DOMAIN = "domain"
    FAMILY = "family"
    MOTIVE = "motive"
    PROTEIN = "protein"
    TURN = "turn"
    SIGNAL = "signal"
    PROPEP = "propeptide"
