import os
from typing import Dict, List, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from pymsaviz import MsaViz
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from .sequence import Sequence
from .standardnumbering import StandardNumbering


class AlignmentResult(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """"""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    consensus: Optional[str] = element(
        description="Consensus sequence of the alignment.",
        default=None,
        tag="consensus",
        json_schema_extra=dict(),
    )

    sequences: List[Sequence] = element(
        description="Sequences of the alignment.",
        default_factory=ListPlus,
        tag="sequences",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    aligned_sequences: List[Sequence] = element(
        description="Aligned sequences as a result of the alignment.",
        default_factory=ListPlus,
        tag="aligned_sequences",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    standard_numbering: Optional[StandardNumbering] = element(
        description="Standard numbering of the aligned sequences.",
        default_factory=StandardNumbering,
        tag="standard_numbering",
        json_schema_extra=dict(),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="72d2203f2e3ce4b319b29fa0d2f146b5eead7b00"
    )

    _raw_xml_data: Dict = PrivateAttr(default_factory=dict)

    @model_validator(mode="after")
    def _parse_raw_xml_data(self):
        for attr, value in self:
            if isinstance(value, (ListPlus, list)) and all(
                isinstance(i, _Element) for i in value
            ):
                self._raw_xml_data[attr] = [elem2dict(i) for i in value]
            elif isinstance(value, _Element):
                self._raw_xml_data[attr] = elem2dict(value)

        return self

    def add_to_sequences(
        self,
        sequence_id: Optional[str] = None,
        sequence: Optional[str] = None,
        id: Optional[str] = None,
        **kwargs,
    ) -> Sequence:
        """
        This method adds an object of type 'Sequence' to attribute sequences

        Args:
            id (str): Unique identifier of the 'Sequence' object. Defaults to 'None'.
            sequence_id (): Identifier of the sequence in the source database. Defaults to None
            sequence (): Molecular sequence.. Defaults to None
        """

        params = {
            "sequence_id": sequence_id,
            "sequence": sequence,
        }

        if id is not None:
            params["id"] = id

        obj = Sequence(**params)

        self.sequences.append(obj)

        return self.sequences[-1]

    def add_to_aligned_sequences(
        self,
        sequence_id: Optional[str] = None,
        sequence: Optional[str] = None,
        id: Optional[str] = None,
        **kwargs,
    ) -> Sequence:
        """
        This method adds an object of type 'Sequence' to attribute aligned_sequences

        Args:
            id (str): Unique identifier of the 'Sequence' object. Defaults to 'None'.
            sequence_id (): Identifier of the sequence in the source database. Defaults to None
            sequence (): Molecular sequence.. Defaults to None
        """

        params = {
            "sequence_id": sequence_id,
            "sequence": sequence,
        }

        if id is not None:
            params["id"] = id

        obj = Sequence(**params)

        self.aligned_sequences.append(obj)

        return self.aligned_sequences[-1]

    def visualize(self):
        """Visualizes the alignment result."""

        # create temp multifasta file
        temp_file = "temp.fasta"
        with open(temp_file, "w") as f:
            for seq in self.aligned_sequences:
                f.write(f">{seq.id}\n{seq.sequence}\n")

        mv = MsaViz(
            temp_file,
            color_scheme="Clustal",
            show_seq_char=False,
            show_count=True,
            sort=True,
            show_consensus=True,
            consensus_color="grey",
        )
        mv.set_plot_params(
            ticks_interval=50, x_unit_size=0.04, show_consensus_char=False
        )
        mv.plotfig()
        os.remove(temp_file)
