from .annotation import Annotation
from .annotationtype import AnnotationType
from .goannotation import GOAnnotation
from .molecule import Molecule
from .protein import Protein
from .pyeedbase import BaseNode
from .reaction import Reaction
from .taxon import Taxon

__all__ = [
    "Annotation",
    "AnnotationType",
    "BaseNode",
    "GOAnnotation",
    "Molecule",
    "Protein",
    "Reaction",
    "Taxon",
]

MODEL_CLASSES: list[type[BaseNode]] = [
    Annotation,
    GOAnnotation,
    Molecule,
    Taxon,
    Protein,
    Reaction,
]
