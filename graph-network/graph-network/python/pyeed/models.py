## This is a generated file. Do not modify it manually!

from __future__ import annotations
from dataclasses import dataclass, field
from dataclasses_json import config, dataclass_json
from typing import List, Optional
from enum import Enum
from uuid import uuid4
from datetime import date, datetime


@dataclass_json
@dataclass
class SequenceRecord:
    sequence: str
    name: Optional[str] = None
    organism: Optional[Organism] = None
    seq_length: Optional[int] = None
    sites: List[Site] = field(default_factory=list)
    regions: List[Region] = field(default_factory=list)
    region_sets: List[RegionSet] = field(default_factory=list)

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:SequenceRecord/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:SequenceRecord",
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
            "sequence": "sio:SIO_000030",
            "id": "sio:SIO_000729",
            "name": "sio:SIO_000116",
            "organism": "sio:SIO_010000",
            "seq_length": "sio:SIO_000041",
        }
    )


    def add_to_sites(
        self,
        positions: List[int]= [],
        **kwargs,
    ):
        params = {
            "positions": positions
        }

        self.sites.append(
            Site(**params)
        )

        return self.sites[-1]


    def add_to_regions(
        self,
        start: Optional[int]= None,
        end: Optional[int]= None,
        **kwargs,
    ):
        params = {
            "start": start,
            "end": end
        }

        self.regions.append(
            Region(**params)
        )

        return self.regions[-1]


    def add_to_region_sets(
        self,
        regions: List[Region]= [],
        **kwargs,
    ):
        params = {
            "regions": regions
        }

        self.region_sets.append(
            RegionSet(**params)
        )

        return self.region_sets[-1]

@dataclass_json
@dataclass
class ProteinRecord:
    structure_id: Optional[str] = None
    coding_sequence: List[Region] = field(default_factory=list)
    ec_number: Optional[str] = None
    mol_weight: Optional[float] = None

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:ProteinRecord/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:ProteinRecord","sio:SIO_010043"
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
            "structure_id": "sio:SIO_000729",
            "coding_sequence": "sio:SIO_001390",
            "ec_number": "edam:data_1011",
            "mol_weight": "edam:data_1505",
        }
    )


    def add_to_coding_sequence(
        self,
        start: Optional[int]= None,
        end: Optional[int]= None,
        **kwargs,
    ):
        params = {
            "start": start,
            "end": end
        }

        self.coding_sequence.append(
            Region(**params)
        )

        return self.coding_sequence[-1]


@dataclass_json
@dataclass
class DNARecord:
    gc_content: Optional[float] = None

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:DNARecord/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:DNARecord","sio:SIO_010008"
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        }
    )


@dataclass_json
@dataclass
class AbstractAnnotation:
    url: Optional[str] = None
    accession_id: Optional[str] = None
    name: Optional[str] = None

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:AbstractAnnotation/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:AbstractAnnotation",
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
            "url": "sio:SIO_000811",
            "accession_id": "sio:SIO_000675",
        }
    )


@dataclass_json
@dataclass
class Site:
    positions: List[int] = field(default_factory=list)

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:Site/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:Site","sio:sio:010049"
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
            "positions": "sio:SIO_000056",
        }
    )


@dataclass_json
@dataclass
class Region:
    start: Optional[int] = None
    end: Optional[int] = None

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:Region/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:Region","sio:SIO_000370"
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
            "start": "sio:SIO_000943",
            "end": "sio:SIO_000953",
        }
    )


@dataclass_json
@dataclass
class RegionSet:
    regions: List[Region] = field(default_factory=list)

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:RegionSet/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:RegionSet","sio:SIO_000370"
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        }
    )


    def add_to_regions(
        self,
        start: Optional[int]= None,
        end: Optional[int]= None,
        **kwargs,
    ):
        params = {
            "start": start,
            "end": end
        }

        self.regions.append(
            Region(**params)
        )

        return self.regions[-1]

@dataclass_json
@dataclass
class Organism:
    taxonomy_id: int
    name: Optional[str] = None
    domain: Optional[str] = None
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    tax_class: Optional[str] = None
    order: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:Organism/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:Organism",
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
            "taxonomy_id": "edam:data_1179",
            "name": "edam:data_2909",
            "kingdom": "edam:data_1044",
            "family": "edam:data_2732",
            "genus": "edam:data_1870",
            "species": "edam:data_1045",
        }
    )


@dataclass_json
@dataclass
class BlastData:
    identity: float = 0.0
    evalue: float = 10.0
    n_hits: int = 100
    substitution_matrix: str = """blosum62"""
    word_size: int = 3
    gap_open: float = 11.0
    gap_extend: float = 1.0
    threshold: float = 11
    db_name: Optional[str] = None

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:BlastData/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:BlastData",
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        }
    )


@dataclass_json
@dataclass
class Sequence:
    sequence_id: Optional[str] = None
    sequence: Optional[str] = None

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:Sequence/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:Sequence",
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        }
    )


@dataclass_json
@dataclass
class AlignmentResult:
    consensus: Optional[str] = None
    sequences: List[Sequence] = field(default_factory=list)
    aligned_sequences: List[Sequence] = field(default_factory=list)
    standard_numbering: Optional[StandardNumbering] = None

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:AlignmentResult/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:AlignmentResult",
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        }
    )


    def add_to_sequences(
        self,
        sequence_id: Optional[str]= None,
        sequence: Optional[str]= None,
        **kwargs,
    ):
        params = {
            "sequence_id": sequence_id,
            "sequence": sequence
        }

        self.sequences.append(
            Sequence(**params)
        )

        return self.sequences[-1]


    def add_to_aligned_sequences(
        self,
        sequence_id: Optional[str]= None,
        sequence: Optional[str]= None,
        **kwargs,
    ):
        params = {
            "sequence_id": sequence_id,
            "sequence": sequence
        }

        self.aligned_sequences.append(
            Sequence(**params)
        )

        return self.aligned_sequences[-1]


@dataclass_json
@dataclass
class PairwiseAlignmentResult:
    score: Optional[float] = None
    identity: Optional[float] = None
    similarity: Optional[float] = None
    gaps: Optional[int] = None
    mismatches: Optional[int] = None

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:PairwiseAlignmentResult/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:PairwiseAlignmentResult",
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        }
    )


@dataclass_json
@dataclass
class StandardNumbering:
    reference_id: Optional[str] = None
    numberd_sequences: List[NumberedSequence] = field(default_factory=list)

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:StandardNumbering/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:StandardNumbering",
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        }
    )


    def add_to_numberd_sequences(
        self,
        numbered_id: Optional[str]= None,
        numbering: List[str]= [],
        **kwargs,
    ):
        params = {
            "numbered_id": numbered_id,
            "numbering": numbering
        }

        self.numberd_sequences.append(
            NumberedSequence(**params)
        )

        return self.numberd_sequences[-1]

@dataclass_json
@dataclass
class NumberedSequence:
    numbered_id: Optional[str] = None
    numbering: List[str] = field(default_factory=list)

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:NumberedSequence/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:NumberedSequence",
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        }
    )


@dataclass_json
@dataclass
class ClustalOmegaResult:
    version: Optional[str] = None

    # JSON-LD fields
    id: str = field(
        metadata=config(field_name="@id"),
        default_factory=lambda: "md:ClustalOmegaResult/" + str(uuid4())
    )
    __type__: list[str] = field(
        metadata=config(field_name="@type"),
        default_factory = lambda: [
            "md:ClustalOmegaResult",
        ],
    )
    __context__: dict[str, str | dict] = field(
        metadata=config(field_name="@context"),
        default_factory = lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        }
    )


class Ontology(Enum):
    ECO = "https://www.evidenceontology.org/term/"
    GO = "https://amigo.geneontology.org/amigo/term/"
    SIO = "http://semanticscience.org/resource/"

class Annotation(Enum):
    ACTIVE_SITE = "http://semanticscience.org/resource/SIO_010041"
    ALLOSTERIC_SITE = "http://semanticscience.org/resource/SIO_010050"
    ALPHAHELIX = "http://semanticscience.org/resource/SIO_010468"
    BETASTRAND = "http://semanticscience.org/resource/SIO_010469"
    BINDING_SITE = "http://semanticscience.org/resource/SIO_010040"
    CODING_SEQ = "http://semanticscience.org/resource/SIO_001276"
    DOMAIN = "http://semanticscience.org/resource/SIO_001379"
    FAMILY = "http://semanticscience.org/resource/SIO_001380"
    MOTIVE = "http://semanticscience.org/resource/SIO_000131"

class SequenceType(Enum):
    DNA = "http://semanticscience.org/resource/SIO_010018"
    PROTEIN = "http://semanticscience.org/resource/SIO_010015"