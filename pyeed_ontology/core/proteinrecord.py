
import validators
from typing import Dict, List, Optional
from pydantic import PrivateAttr, field_validator, model_validator
from uuid import uuid4
from pydantic_xml import attr, element
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from .region import Region
from .position import Position
from .doi import DOI
from .sequencerecord import SequenceRecord


@forge_signature
class ProteinRecord(SequenceRecord, search_mode="unordered"):
    """A protein sequence and associated metadata."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    accession_id: Optional[str] = element(
        description="Accession ID of the protein sequence.",
        default=None,
        tag="accession_id",
        json_schema_extra=dict(term="http://edamontology.org/data_2091"),
    )

    name: Optional[str] = element(
        description="Name of the protein sequence.",
        default=None,
        tag="name",
        json_schema_extra=dict(term="http://edamontology.org/data_1009"),
    )

    sequence: str = element(
        description="Amino acid sequence of the protein.",
        tag="sequence",
        json_schema_extra=dict(term="http://edamontology.org/data_2976"),
    )

    families: List[Region] = element(
        description="Family of the protein",
        default_factory=ListPlus,
        tag="families",
        json_schema_extra=dict(multiple=True, term="http://edamontology.org/data_1131"),
    )

    domains: List[Region] = element(
        description="Domains of the protein",
        default_factory=ListPlus,
        tag="domains",
        json_schema_extra=dict(multiple=True, term="http://edamontology.org/data_1468"),
    )

    sites: List[Position] = element(
        description="Annotations of different sites",
        default_factory=ListPlus,
        tag="sites",
        json_schema_extra=dict(multiple=True),
    )

    coding_sequence: List[Region] = element(
        description="Defines the coding sequence of the protein",
        default_factory=ListPlus,
        tag="coding_sequence",
        json_schema_extra=dict(multiple=True),
    )

    ec_number: Optional[str] = element(
        description="An Enzyme Commission (EC) number of an enzyme.",
        default=None,
        tag="ec_number",
        json_schema_extra=dict(term="http://edamontology.org/data_1011"),
    )

    mol_weight: Optional[float] = element(
        description="Calculated molecular weight of the protein",
        default=None,
        tag="mol_weight",
        json_schema_extra=dict(),
    )

    pdb_id: Optional[str] = element(
        description="Protein Data Bank (PDB) identifier.",
        default=None,
        tag="pdb_id",
        json_schema_extra=dict(),
    )

    alphafold_id: Optional[str] = element(
        description="AlphaFold identifier.",
        default=None,
        tag="alphafold_id",
        json_schema_extra=dict(),
    )

    annotations_: List[str] = element(
        tag="annotations_",
        alias="@type",
        description="Annotation of the given object.",
        default=["ProteinRecord"],
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

    @field_validator("annotations_")
    @classmethod
    def _validate_annotation(cls, annotations):
        """Check if the annotation that has been set is a valid URL."""
        for annotation in annotations:
            if not validators.url(annotation) and annotation != cls.__name__:
                raise ValueError(f"Invalid Annotation URL: {annotation}")
        return annotations

    def add_to_families(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        accession_id: Optional[str] = None,
        name: Optional[str] = None,
        publications: List[DOI] = ListPlus(),
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Region:
        """
        This method adds an object of type 'Region' to attribute families

        Args:
            id (str): Unique identifier of the 'Region' object. Defaults to 'None'.
            start (): Start position in the amino acid sequence of the annotated region.. Defaults to None
            end (): End position position in the amino acid sequence of the annotated region.. Defaults to None
            accession_id (): Accession ID of the annotation.. Defaults to None
            name (): A name of a sequence feature, e.g. the name of a feature. Defaults to None
            publications (): DOI of the publication(s) referenced in the annotation.. Defaults to ListPlus()
        """
        params = {
            "start": start,
            "end": end,
            "accession_id": accession_id,
            "name": name,
            "publications": publications,
        }
        if id is not None:
            params["id"] = id
        obj = Region(**params)
        if annotation:
            obj.annotations_.append(annotation)
        self.families.append(obj)
        return self.families[-1]

    def add_to_domains(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        accession_id: Optional[str] = None,
        name: Optional[str] = None,
        publications: List[DOI] = ListPlus(),
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Region:
        """
        This method adds an object of type 'Region' to attribute domains

        Args:
            id (str): Unique identifier of the 'Region' object. Defaults to 'None'.
            start (): Start position in the amino acid sequence of the annotated region.. Defaults to None
            end (): End position position in the amino acid sequence of the annotated region.. Defaults to None
            accession_id (): Accession ID of the annotation.. Defaults to None
            name (): A name of a sequence feature, e.g. the name of a feature. Defaults to None
            publications (): DOI of the publication(s) referenced in the annotation.. Defaults to ListPlus()
        """
        params = {
            "start": start,
            "end": end,
            "accession_id": accession_id,
            "name": name,
            "publications": publications,
        }
        if id is not None:
            params["id"] = id
        obj = Region(**params)
        if annotation:
            obj.annotations_.append(annotation)
        self.domains.append(obj)
        return self.domains[-1]

    def add_to_sites(
        self,
        position: Optional[int] = None,
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Position:
        """
        This method adds an object of type 'Position' to attribute sites

        Args:
            id (str): Unique identifier of the 'Position' object. Defaults to 'None'.
            position (): Position in a sequence.. Defaults to None
        """
        params = {"position": position}
        if id is not None:
            params["id"] = id
        obj = Position(**params)
        if annotation:
            obj.annotations_.append(annotation)
        self.sites.append(obj)
        return self.sites[-1]

    def add_to_coding_sequence(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        accession_id: Optional[str] = None,
        name: Optional[str] = None,
        publications: List[DOI] = ListPlus(),
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Region:
        """
        This method adds an object of type 'Region' to attribute coding_sequence

        Args:
            id (str): Unique identifier of the 'Region' object. Defaults to 'None'.
            start (): Start position in the amino acid sequence of the annotated region.. Defaults to None
            end (): End position position in the amino acid sequence of the annotated region.. Defaults to None
            accession_id (): Accession ID of the annotation.. Defaults to None
            name (): A name of a sequence feature, e.g. the name of a feature. Defaults to None
            publications (): DOI of the publication(s) referenced in the annotation.. Defaults to ListPlus()
        """
        params = {
            "start": start,
            "end": end,
            "accession_id": accession_id,
            "name": name,
            "publications": publications,
        }
        if id is not None:
            params["id"] = id
        obj = Region(**params)
        if annotation:
            obj.annotations_.append(annotation)
        self.coding_sequence.append(obj)
        return self.coding_sequence[-1]

    @classmethod
    def _validate_annotation(cls, annotations):
        """Check if the annotation that has been set is a valid URL."""

        for annotation in annotations:
            if not validators.url(annotation) and annotation != cls.__name__:
                raise ValueError(f"Invalid Annotation URL: {annotation}")

        return annotations

    def add_to_families(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        accession_id: Optional[str] = None,
        name: Optional[str] = None,
        publications: List[DOI] = ListPlus(),
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Region:
        """
        This method adds an object of type 'Region' to attribute families

        Args:
            id (str): Unique identifier of the 'Region' object. Defaults to 'None'.
            start (): Start position in the amino acid sequence of the annotated region.. Defaults to None
            end (): End position position in the amino acid sequence of the annotated region.. Defaults to None
            accession_id (): Accession ID of the annotation.. Defaults to None
            name (): A name of a sequence feature, e.g. the name of a feature. Defaults to None
            publications (): DOI of the publication(s) referenced in the annotation.. Defaults to ListPlus()
        """

        params = {
            "start": start,
            "end": end,
            "accession_id": accession_id,
            "name": name,
            "publications": publications,
        }

        if id is not None:
            params["id"] = id

        obj = Region(**params)

        if annotation:
            obj.annotations_.append(annotation)

        self.families.append(obj)

        return self.families[-1]

    def add_to_domains(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        accession_id: Optional[str] = None,
        name: Optional[str] = None,
        publications: List[DOI] = ListPlus(),
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Region:
        """
        This method adds an object of type 'Region' to attribute domains

        Args:
            id (str): Unique identifier of the 'Region' object. Defaults to 'None'.
            start (): Start position in the amino acid sequence of the annotated region.. Defaults to None
            end (): End position position in the amino acid sequence of the annotated region.. Defaults to None
            accession_id (): Accession ID of the annotation.. Defaults to None
            name (): A name of a sequence feature, e.g. the name of a feature. Defaults to None
            publications (): DOI of the publication(s) referenced in the annotation.. Defaults to ListPlus()
        """

        params = {
            "start": start,
            "end": end,
            "accession_id": accession_id,
            "name": name,
            "publications": publications,
        }

        if id is not None:
            params["id"] = id

        obj = Region(**params)

        if annotation:
            obj.annotations_.append(annotation)

        self.domains.append(obj)

        return self.domains[-1]

    def add_to_sites(
        self,
        position: Optional[int] = None,
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Position:
        """
        This method adds an object of type 'Position' to attribute sites

        Args:
            id (str): Unique identifier of the 'Position' object. Defaults to 'None'.
            position (): Position in a sequence.. Defaults to None
        """

        params = {
            "position": position,
        }

        if id is not None:
            params["id"] = id

        obj = Position(**params)

        if annotation:
            obj.annotations_.append(annotation)

        self.sites.append(obj)

        return self.sites[-1]

    def add_to_coding_sequence(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        accession_id: Optional[str] = None,
        name: Optional[str] = None,
        publications: List[DOI] = ListPlus(),
        id: Optional[str] = None,
        annotation: Optional[str] = None,
        **kwargs,
    ) -> Region:
        """
        This method adds an object of type 'Region' to attribute coding_sequence

        Args:
            id (str): Unique identifier of the 'Region' object. Defaults to 'None'.
            start (): Start position in the amino acid sequence of the annotated region.. Defaults to None
            end (): End position position in the amino acid sequence of the annotated region.. Defaults to None
            accession_id (): Accession ID of the annotation.. Defaults to None
            name (): A name of a sequence feature, e.g. the name of a feature. Defaults to None
            publications (): DOI of the publication(s) referenced in the annotation.. Defaults to ListPlus()
        """

        params = {
            "start": start,
            "end": end,
            "accession_id": accession_id,
            "name": name,
            "publications": publications,
        }

        if id is not None:
            params["id"] = id

        obj = Region(**params)

        if annotation:
            obj.annotations_.append(annotation)

        self.coding_sequence.append(obj)

        return self.coding_sequence[-1]
