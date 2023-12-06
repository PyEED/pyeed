import sdRDM

from typing import Optional, Union, List
from pydantic import PrivateAttr, Field, validator
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .proteinsitetype import ProteinSiteType
from .dnaregiontype import DNARegionType
from .proteinregiontype import ProteinRegionType
from .proteinregion import ProteinRegion
from .dnaregion import DNARegion
from .site import Site
from .organism import Organism
from .substrate import Substrate
from .span import Span
from .author import Author
from .proteininfo import ProteinInfo
from .citation import Citation


@forge_signature
class PairwiseAlignment(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("pairwisealignmentINDEX"),
        xml="@id",
    )

    reference_seq: Optional[ProteinInfo] = Field(
        default=None,
        description="Protein sequence used as reference",
        alias="reference",
    )

    query_seq: Optional[ProteinInfo] = Field(
        default=None,
        description="Protein sequence used as query",
        alias="query",
    )

    score: Optional[float] = Field(
        default=None,
        description="Alignment score",
    )

    identity: Optional[float] = Field(
        default=None,
        description="Ration of identical residues in the alignment",
    )

    similarity: Optional[float] = Field(
        default=None,
        description="Ration of similar residues in the alignment",
    )

    gaps: Optional[int] = Field(
        default=None,
        description="Number of gaps in the alignment",
    )

    mismatches: Optional[int] = Field(
        default=None,
        description="Number of mismatches in the alignment",
    )
