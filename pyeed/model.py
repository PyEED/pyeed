## This is a generated file. Do not modify it manually!

from __future__ import annotations

from enum import Enum
from typing import Generic, Optional, TypeVar
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field

# Filter Wrapper definition used to filter a list of objects
# based on their attributes
Cls = TypeVar("Cls")


class FilterWrapper(Generic[Cls]):
    """Wrapper class to filter a list of objects based on their attributes"""

    def __init__(self, collection: list[Cls], **kwargs):
        self.collection = collection
        self.kwargs = kwargs

    def filter(self) -> list[Cls]:
        for key, value in self.kwargs.items():
            self.collection = [
                item for item in self.collection if self._fetch_attr(key, item) == value
            ]
        return self.collection

    def _fetch_attr(self, name: str, item: Cls):
        try:
            return getattr(item, name)
        except AttributeError:
            raise AttributeError(f"{item} does not have attribute {name}")


# JSON-LD Helper Functions
def add_namespace(obj, prefix: str | None, iri: str | None):
    """Adds a namespace to the JSON-LD context

    Args:
        prefix (str): The prefix to add
        iri (str): The IRI to add
    """
    if prefix is None and iri is None:
        return
    elif prefix and iri is None:
        raise ValueError("If prefix is provided, iri must also be provided")
    elif iri and prefix is None:
        raise ValueError("If iri is provided, prefix must also be provided")

    obj.ld_context[prefix] = iri  # type: ignore


def validate_prefix(term: str | dict, prefix: str):
    """Validates that a term is prefixed with a given prefix

    Args:
        term (str): The term to validate
        prefix (str): The prefix to validate against

    Returns:
        bool: True if the term is prefixed with the prefix, False otherwise
    """

    if isinstance(term, dict) and not term["@id"].startswith(prefix + ":"):
        raise ValueError(f"Term {term} is not prefixed with {prefix}")
    elif isinstance(term, str) and not term.startswith(prefix + ":"):
        raise ValueError(f"Term {term} is not prefixed with {prefix}")


# Model Definitions


class ProteinRecord(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
        use_enum_values=True,
    )  # type: ignore

    id: str
    sequence: str
    name: Optional[str] = Field(default=None)
    organism: Optional[Organism] = Field(default=None)
    embedding: list[float] = Field(default_factory=list)
    seq_length: Optional[int] = Field(default=None)
    nucleotide_id: Optional[str] = Field(default=None)
    locus_tag: Optional[str] = Field(default=None)
    sites: list[Site] = Field(default_factory=list)
    regions: list[Region] = Field(default_factory=list)
    structure_ids: list[str] = Field(default_factory=list)
    ec_number: Optional[str] = Field(default=None)
    mol_weight: Optional[float] = Field(default=None)
    annotations: list[Annotation] = Field(default_factory=list)
    go_terms: list[str] = Field(default_factory=list)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:ProteinRecord/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: ["md:ProteinRecord", "sio:SIO_010043"],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
            "id": "sio:SIO_000729",
            "sequence": "sio:SIO_000030",
            "name": "sio:SIO_000116",
            "organism": "sio:SIO_010000",
            "seq_length": "sio:SIO_000041",
            "structure_ids": "sio:SIO_000729",
            "ec_number": "edam:data_1011",
            "mol_weight": "edam:data_1505",
        },
    )

    def filter_sites(self, **kwargs) -> list[Site]:
        """Filters the sites attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[Site]: The filtered list of Site objects
        """

        return FilterWrapper[Site](self.sites, **kwargs).filter()

    def filter_regions(self, **kwargs) -> list[Region]:
        """Filters the regions attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[Region]: The filtered list of Region objects
        """

        return FilterWrapper[Region](self.regions, **kwargs).filter()

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)

    def add_to_sites(
        self,
        annotation: Annotation,
        name: Optional[str] = None,
        positions: list[int] = [],
        **kwargs,
    ):
        params = {"annotation": annotation, "name": name, "positions": positions}

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.sites.append(Site(**params))

        return self.sites[-1]

    def add_to_regions(
        self,
        id: str,
        annotation: Annotation,
        start: int,
        end: int,
        name: Optional[str] = None,
        **kwargs,
    ):
        params = {
            "id": id,
            "annotation": annotation,
            "start": start,
            "end": end,
            "name": name,
        }

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.regions.append(Region(**params))

        return self.regions[-1]


class DNARecord(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
        use_enum_values=True,
    )  # type: ignore

    id: str
    sequence: str
    name: Optional[str] = Field(default=None)
    organism: Optional[Organism] = Field(default=None)
    seq_length: Optional[int] = Field(default=None)
    embedding: list[float] = Field(default_factory=list)
    sites: list[Site] = Field(default_factory=list)
    regions: list[Region] = Field(default_factory=list)
    gc_content: Optional[float] = Field(default=None)
    annotations: list[Annotation] = Field(default_factory=list)
    go_terms: list[str] = Field(default_factory=list)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:DNARecord/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: ["md:DNARecord", "sio:SIO_010008"],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
            "id": "sio:SIO_000729",
            "sequence": "sio:SIO_000030",
            "name": "sio:SIO_000116",
            "organism": "sio:SIO_010000",
            "seq_length": "sio:SIO_000041",
        },
    )

    def filter_sites(self, **kwargs) -> list[Site]:
        """Filters the sites attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[Site]: The filtered list of Site objects
        """

        return FilterWrapper[Site](self.sites, **kwargs).filter()

    def filter_regions(self, **kwargs) -> list[Region]:
        """Filters the regions attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[Region]: The filtered list of Region objects
        """

        return FilterWrapper[Region](self.regions, **kwargs).filter()

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)

    def add_to_sites(
        self,
        annotation: Annotation,
        name: Optional[str] = None,
        positions: list[int] = [],
        **kwargs,
    ):
        params = {"annotation": annotation, "name": name, "positions": positions}

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.sites.append(Site(**params))

        return self.sites[-1]

    def add_to_regions(
        self,
        id: str,
        annotation: Annotation,
        start: int,
        end: int,
        name: Optional[str] = None,
        **kwargs,
    ):
        params = {
            "id": id,
            "annotation": annotation,
            "start": start,
            "end": end,
            "name": name,
        }

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.regions.append(Region(**params))

        return self.regions[-1]


class Site(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
        use_enum_values=True,
    )  # type: ignore

    annotation: Annotation
    name: Optional[str] = Field(default=None)
    positions: list[int] = Field(default_factory=list)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id", default_factory=lambda: "md:Site/" + str(uuid4())
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: ["md:Site", "sio:sio:010049"],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
            "positions": "sio:SIO_000056",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class Region(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
        use_enum_values=True,
    )  # type: ignore

    id: str
    annotation: Annotation
    start: int
    end: int
    name: Optional[str] = Field(default=None)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id", default_factory=lambda: "md:Region/" + str(uuid4())
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:Region",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
            "start": "sio:SIO_000943",
            "end": "sio:SIO_000953",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class Organism(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
        use_enum_values=True,
    )  # type: ignore

    taxonomy_id: int
    name: Optional[str] = Field(default=None)
    domain: Optional[str] = Field(default=None)
    kingdom: Optional[str] = Field(default=None)
    phylum: Optional[str] = Field(default=None)
    tax_class: Optional[str] = Field(default=None)
    order: Optional[str] = Field(default=None)
    family: Optional[str] = Field(default=None)
    genus: Optional[str] = Field(default=None)
    species: Optional[str] = Field(default=None)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id", default_factory=lambda: "md:Organism/" + str(uuid4())
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:Organism",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
            "taxonomy_id": "edam:data_1179",
            "name": "edam:data_2909",
            "kingdom": "edam:data_1044",
            "family": "edam:data_2732",
            "genus": "edam:data_1870",
            "species": "edam:data_1045",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class BlastData(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
        use_enum_values=True,
    )  # type: ignore

    identity: float = 0.0
    evalue: float = 10.0
    n_hits: int = 100
    substitution_matrix: str = """blosum62"""
    word_size: int = 3
    gap_open: float = 11.0
    gap_extend: float = 1.0
    threshold: float = 11
    db_name: Optional[str] = Field(default=None)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:BlastData/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:BlastData",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class Sequence(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
        use_enum_values=True,
    )  # type: ignore

    sequence_id: Optional[str] = Field(default=None)
    sequence: Optional[str] = Field(default=None)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id", default_factory=lambda: "md:Sequence/" + str(uuid4())
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:Sequence",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class AlignmentResult(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
        use_enum_values=True,
    )  # type: ignore

    consensus: Optional[str] = Field(default=None)
    sequences: list[Sequence] = Field(default_factory=list)
    aligned_sequences: list[Sequence] = Field(default_factory=list)
    standard_numbering: Optional[StandardNumbering] = Field(default=None)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:AlignmentResult/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:AlignmentResult",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        },
    )

    def filter_sequences(self, **kwargs) -> list[Sequence]:
        """Filters the sequences attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[Sequence]: The filtered list of Sequence objects
        """

        return FilterWrapper[Sequence](self.sequences, **kwargs).filter()

    def filter_aligned_sequences(self, **kwargs) -> list[Sequence]:
        """Filters the aligned_sequences attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[Sequence]: The filtered list of Sequence objects
        """

        return FilterWrapper[Sequence](self.aligned_sequences, **kwargs).filter()

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)

    def add_to_sequences(
        self,
        sequence_id: Optional[str] = None,
        sequence: Optional[str] = None,
        **kwargs,
    ):
        params = {"sequence_id": sequence_id, "sequence": sequence}

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.sequences.append(Sequence(**params))

        return self.sequences[-1]

    def add_to_aligned_sequences(
        self,
        sequence_id: Optional[str] = None,
        sequence: Optional[str] = None,
        **kwargs,
    ):
        params = {"sequence_id": sequence_id, "sequence": sequence}

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.aligned_sequences.append(Sequence(**params))

        return self.aligned_sequences[-1]


class PairwiseAlignmentResult(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
        use_enum_values=True,
    )  # type: ignore

    score: Optional[float] = Field(default=None)
    identity: Optional[float] = Field(default=None)
    similarity: Optional[float] = Field(default=None)
    gaps: Optional[int] = Field(default=None)
    mismatches: Optional[int] = Field(default=None)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:PairwiseAlignmentResult/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:PairwiseAlignmentResult",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class StandardNumbering(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
        use_enum_values=True,
    )  # type: ignore

    reference_id: Optional[str] = Field(default=None)
    numberd_sequences: list[NumberedSequence] = Field(default_factory=list)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:StandardNumbering/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:StandardNumbering",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        },
    )

    def filter_numberd_sequences(self, **kwargs) -> list[NumberedSequence]:
        """Filters the numberd_sequences attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[NumberedSequence]: The filtered list of NumberedSequence objects
        """

        return FilterWrapper[NumberedSequence](
            self.numberd_sequences, **kwargs
        ).filter()

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)

    def add_to_numberd_sequences(
        self,
        numbered_id: Optional[str] = None,
        numbering: list[str] = [],
        **kwargs,
    ):
        params = {"numbered_id": numbered_id, "numbering": numbering}

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.numberd_sequences.append(NumberedSequence(**params))

        return self.numberd_sequences[-1]


class NumberedSequence(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
        use_enum_values=True,
    )  # type: ignore

    numbered_id: Optional[str] = Field(default=None)
    numbering: list[str] = Field(default_factory=list)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:NumberedSequence/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:NumberedSequence",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "http://mdmodel.net/",
            "sio": "http://semanticscience.org/resource/",
            "edam": "http://edamontology.org/",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class Annotation(Enum):
    ACTIVE_SITE = "http://semanticscience.org/resource/SIO_010041"
    ALLOSTERIC_SITE = "http://semanticscience.org/resource/SIO_010050"
    ALPHAHELIX = "http://semanticscience.org/resource/SIO_010468"
    BETASTRAND = "http://semanticscience.org/resource/SIO_010469"
    BINDING_SITE = "http://semanticscience.org/resource/SIO_010040"
    CODING_SEQ = "http://semanticscience.org/resource/SIO_001276"
    DNA = "http://semanticscience.org/resource/SIO_010018"
    DOMAIN = "http://semanticscience.org/resource/SIO_001379"
    FAMILY = "http://semanticscience.org/resource/SIO_001380"
    MOTIVE = "http://semanticscience.org/resource/SIO_000131"
    PROTEIN = "http://semanticscience.org/resource/SIO_010015"
