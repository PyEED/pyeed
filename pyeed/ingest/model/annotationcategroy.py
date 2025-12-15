from enum import Enum


class AnnotationCategory(str, Enum):
    """Protein and DNA annotation types."""

    SITE = "site"
    ACTIVE_SITE = "active_site"
    BINDING_SITE = "binding_site"
    DOMAIN = "domain"
    FAMILY = "family"
    SUPERFAMILY = "superfamily"
