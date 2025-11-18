from .annotation import Annotation
from .annotationtype import AnnotationType
from .goannotation import GOAnnotation
from .molecule import Molecule
from .protein import Protein
from .pyeedbase import PyeedBase
from .reaction import Reaction
from .taxon import Taxon

__all__ = [
    "Annotation",
    "AnnotationType",
    "Embedding",
    "GOAnnotation",
    "Molecule",
    "Protein",
    "PyeedBase",
    "Reaction",
    "Taxon",
]

MODEL_CLASSES: list[type[PyeedBase]] = [
    Annotation,
    GOAnnotation,
    Molecule,
    Taxon,
    Protein,
    Reaction,
]
