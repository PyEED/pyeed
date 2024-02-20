import sdRDM

from typing import Dict, List, Optional
from pydantic import PrivateAttr, model_validator
from uuid import uuid4
from pydantic_xml import attr, element, wrapped
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from .dnaregion import DNARegion
from .dnaregiontype import DNARegionType
from .span import Span
from ..ncbi.seq_io import get_ncbi_entry, _seqio_to_dna_info


@forge_signature
class DNAInfo(
    sdRDM.DataModel,
    nsmap={
        "": "https://github.com/PyEED/pyeed@3b002efa6bf51e951767d8a7749ebad563897cb8#DNAInfo"
    },
):
    """Description of a nucleotide sequence 🧬"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    regions: List[DNARegion] = wrapped(
        "regions",
        element(
            description=(
                "Defines regions within the nucleotide sequence that code for the"
                " protein sequence"
            ),
            default_factory=ListPlus,
            tag="DNARegion",
            json_schema_extra=dict(multiple=True),
        ),
    )
    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="3b002efa6bf51e951767d8a7749ebad563897cb8"
    )
    _raw_xml_data: Dict = PrivateAttr(default_factory=dict)

    @model_validator(mode="after")
    def _parse_raw_xml_data(self):
        for attr, value in self:
            if isinstance(value, (ListPlus, list)) and all(
                (isinstance(i, _Element) for i in value)
            ):
                self._raw_xml_data[attr] = [elem2dict(i) for i in value]
            elif isinstance(value, _Element):
                self._raw_xml_data[attr] = elem2dict(value)
        return self

    def add_to_regions(
        self,
        type: Optional[DNARegionType] = None,
        name: Optional[str] = None,
        spans: List[Span] = ListPlus(),
        note: Optional[str] = None,
        cross_reference: Optional[str] = None,
        id: Optional[str] = None,
    ) -> DNARegion:
        """
        This method adds an object of type 'DNARegion' to attribute regions

        Args:
            id (str): Unique identifier of the 'DNARegion' object. Defaults to 'None'.
            type (): Type of the region within the nucleotide sequence. Defaults to None
            name (): Name of the annotation. Defaults to None
            spans (): Spans of the region. E.g. multiple exons of a gene. Defaults to ListPlus()
            note (): Information found in 'note' of an ncbi entry. Defaults to None
            cross_reference (): Database cross reference. Defaults to None
        """
        params = {
            "type": type,
            "name": name,
            "spans": spans,
            "note": note,
            "cross_reference": cross_reference,
        }
        if id is not None:
            params["id"] = id
        self.regions.append(DNARegion(**params))
        return self.regions[-1]

    @classmethod
    def from_ncbi(cls, accession_id: str) -> "DNAInfo":
        seq_record = get_ncbi_entry(accession_id=accession_id, database="nucleotide")
        return _seqio_to_dna_info(cls, seq_record)

    # def extract_nucleotide_seq(self, dna_region: DNARegion):
    #     """Handel nucleotide SeqIO entry and map it to `NucleotideSequence`"""

    #     for feature in entry.features:
    #         if feature.type == "CDS":
    #             if isinstance(feature.location, CompoundLocation):
    #                 locations = set()
    #                 parts = feature.location.parts
    #                 for part in parts:
    #                     # TODO: investigate reason for +1
    #                     locations.add(int(part.start) + 1)
    #                     locations.add(int(part.end))

    #             if isinstance(feature.location, FeatureLocation):
    #                 locations = {int(feature.location.start) + 1, int(feature.location.end)}

    #             if feature_regions == locations:
    #                 nucleotide_sequence.sequence = str(feature.location.extract(entry.seq))
