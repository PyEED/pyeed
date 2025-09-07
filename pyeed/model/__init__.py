from .annotation import Annotation
from .annotationtype import AnnotationType
from .embedding import Embedding
from .goannotation import GOAnnotation
from .molecule import Molecule
from .organism import Organism
from .protein import Protein
from .pyeedbase import PyeedBase
from .reaction import Reaction

__all__ = [
    "Annotation",
    "AnnotationType",
    "Embedding",
    "GOAnnotation",
    "Molecule",
    "Organism",
    "Protein",
    "PyeedBase",
    "Reaction",
]

MODEL_CLASSES: list[type[PyeedBase]] = [
    Annotation,
    Embedding,
    GOAnnotation,
    Molecule,
    Organism,
    Protein,
    Reaction,
]
