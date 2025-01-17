import asyncio
import warnings
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional
from uuid import uuid4

import nest_asyncio
from Bio.Blast import NCBIXML
from IPython.display import clear_output
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from pyeed.fetch.blast import Blast, BlastProgram, NCBIDataBase
from pyeed.fetch.proteinfetcher import ProteinFetcher
from rich.console import Console
from rich.status import Status
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from .dnarecord import DNARecord
from .region import Region
from .sequencerecord import SequenceRecord


class ProteinRecord(
    SequenceRecord,
):
    """A protein sequence and associated metadata."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    structure_id: Optional[str] = element(
        description="Protein Data Bank (PDB) identifier.",
        default=None,
        tag="structure_id",
        json_schema_extra=dict(
            term="http://semanticscience.org/resource/SIO_000729",
        ),
    )

    coding_sequence: List[Region] = element(
        description="Defines the coding sequence of the protein",
        default_factory=ListPlus,
        tag="coding_sequence",
        json_schema_extra=dict(
            multiple=True,
            term="http://semanticscience.org/resource/SIO_001390",
        ),
    )

    ec_number: Optional[str] = element(
        description="An Enzyme Commission (EC) number of an enzyme.",
        default=None,
        tag="ec_number",
        json_schema_extra=dict(
            term="http://edamontology.org/data_1011",
        ),
    )

    mol_weight: Optional[float] = element(
        description="Calculated molecular weight of the protein based on the sequence.",
        default=None,
        tag="mol_weight",
        json_schema_extra=dict(
            term="http://edamontology.org/data_1505",
        ),
    )

    _repo: Optional[str] = PrivateAttr(default="https://github.com/PyEED/pyeed")
    _commit: Optional[str] = PrivateAttr(
        default="b926bfec3aa1ec45a5614cf6ac4a546252dd384c"
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

    def add_to_coding_sequence(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        url: Optional[str] = None,
        accession_id: Optional[str] = None,
        name: Optional[str] = None,
        id: Optional[str] = None,
        **kwargs,
    ) -> Region:
        """
        This method adds an object of type 'Region' to attribute coding_sequence

        Args:
            id (str): Unique identifier of the 'Region' object. Defaults to 'None'.
            start (): Start position of the site.. Defaults to None
            end (): End position of the site.. Defaults to None
            url (): URI of the annotation.. Defaults to None
            accession_id (): Accession ID of the annotation.. Defaults to None
            name (): A name of a sequence feature, e.g. the name of a feature. Defaults to None
        """

        params = {
            "start": start,
            "end": end,
            "url": url,
            "accession_id": accession_id,
            "name": name,
        }

        if id is not None:
            params["id"] = id

        obj = Region(**params)

        self.coding_sequence.append(obj)

        return self.coding_sequence[-1]

    @classmethod
    def get_id(cls, protein_id: str) -> "ProteinRecord":
        """
        This method creates a 'ProteinRecord' object from a given protein accession ID.

        Args:
            protein_id (str): ID of the protein in NCBI or UniProt database.

        Returns:
            ProteinRecord: 'ProteinRecord' with information of the corresponding protein_id.
        """

        import nest_asyncio

        nest_asyncio.apply()

        if isinstance(protein_id, list) and all(isinstance(x, str) for x in protein_id):
            warnings.warn("For getting multiple sequences by ID use `get_ids` instead.")
            return cls.get_ids(protein_id)

        sequences = asyncio.run(ProteinFetcher(ids=[protein_id]).fetch(quiet=False))[0]
        clear_output()
        return sequences

    @classmethod
    def get_ids(cls, accession_ids: List[str]) -> List["ProteinRecord"]:
        """Creates a list of 'ProteinRecord' objects from a list of protein accession IDs.

        Returns:
            List[ProteinRecord]: A list of 'ProteinRecord' objects representing the protein sequences found.
        """

        import nest_asyncio

        nest_asyncio.apply()

        return asyncio.run(
            ProteinFetcher(ids=accession_ids).fetch(force_terminal=False)
        )

    @classmethod
    def from_sequence(
        cls,
        sequence: str,
        exact_match: bool = True,
        database: str = "nr",
        matrix: str = "BLOSUM62",
    ):
        """
        Creates a 'ProteinInfo' object from a given protein sequence by
        performing a BLAST search on NCBI server.

        Args:
            sequence (str): The protein sequence to search for.
            exact_match (bool, optional): If True, only exact matches will be considered.
                If False, approximate matches will also be included. Defaults to True.
            database (str, optional): The database to search against. Must be one of
                the supported databases: 'nr', 'swissprot', 'pdb', 'refseq_protein'.
                Defaults to 'nr'.

        Returns:
            ProteinInfo: A 'ProteinInfo' object representing the protein sequence
                found in the database.

        Raises:
            AssertionError: If the specified database is not supported.
        """

        from pyeed.fetch.blast import BlastProgram

        nest_asyncio.apply()

        assert (
            database in NCBIDataBase
        ), f"Database needs to be one of {NCBIDataBase.__members__.keys()}"

        identity = 1 if exact_match else 0

        blaster = Blast(
            query=sequence,
            n_hits=1,
            identity=identity,
            matrix=matrix,
        )

        with Status("Running BLAST", console=Console(force_terminal=False)) as status:
            result = asyncio.run(
                blaster.async_run(
                    NCBIDataBase.NR.value,
                    BlastProgram.BLASTP.value,
                )
            )
            clear_output()

            accession = blaster.extract_accession(result)

            status.update("Fetching protein data")

            if accession:
                return asyncio.run(
                    ProteinFetcher(ids=accession).fetch(force_terminal=False)
                )[0]

        return

    def ncbi_blast(
        self,
        n_hits: int,
        e_value: float = 10.0,
        db: str = "swissprot",
        matrix: str = "BLOSUM62",
        identity: float = 0.0,
        **kwargs,
    ) -> List["ProteinRecord"]:
        """
        Runs a BLAST search using the NCBI BLAST service to find similar protein sequences.

        Args:
            n_hits (int): The number of hits to retrieve.
            e_value (float, optional): The maximum E-value threshold for reporting hits. Defaults to 10.0.
            db (str, optional): The database to search against. Defaults to "swissprot".
            matrix (str, optional): The substitution matrix to use. Defaults to "BLOSUM62".
            identity (float, optional): The minimum sequence identity threshold for reporting hits. Defaults to 0.0.
            **kwargs: Additional keyword arguments.

        Returns:
            List[ProteinInfo]: A list of ProteinInfo objects representing the similar protein sequences found.

        Raises:
            AssertionError: If the specified database is not supported.

        Example:
            protein_info = ProteinInfo()
            similar_proteins = protein_info.ncbi_blast(n_hits=10, e_value=0.001, database="swissprot")
        """

        import nest_asyncio
        from pyeed.fetch.blast import NCBIDataBase

        nest_asyncio.apply()

        assert (
            db in NCBIDataBase
        ), f"Database needs to be one of {NCBIDataBase.__members__.keys()}"

        program = BlastProgram.BLASTP.value
        executor = ThreadPoolExecutor(max_workers=1)
        blaster = Blast(
            query=self.sequence,
            n_hits=n_hits,
            evalue=e_value,
            matrix=matrix,
            identity=identity,
        )

        with Status(
            "Running BLAST", console=Console(force_terminal=False, force_jupyter=True)
        ):
            result = asyncio.run(blaster.async_run(db, program, executor))
            clear_output()

        accessions = blaster.extract_accession(result)

        return self.get_ids(accessions)

    # def blastp(
    #     self,
    #     db_path: str,
    #     identity: float = 0,
    #     evalue: float = 10,
    #     n_hits: int = 500,
    #     subst_matrix: str = "BLOSUM62",
    #     word_size: int = 3,
    #     gapopen: int = 11,
    #     gapextend: int = 1,
    #     threshold: int = 11,
    #     n_cores: int = os.cpu_count(),
    #     ncbi_key: str = None,
    #     email: str = None,
    # ):
    #     blaster = Blastp(
    #         _db_path=db_path,
    #         identity=identity,
    #         evalue=evalue,
    #         n_hits=n_hits,
    #         subst_matrix=subst_matrix,
    #         word_size=word_size,
    #         gapopen=gapopen,
    #         gapextend=gapextend,
    #         threshold=threshold,
    #         n_cores=n_cores,
    #         ncbi_key=ncbi_key,
    #     )

    #     command = blaster.setup_command()
    #     accession_ids = blaster.run_container(
    #         command=command, data=self._fasta_string()
    #     )
    #     protein_infos = ProteinRecord.get_ids(accession_ids, email, ncbi_key)
    #     protein_infos.insert(0, self)
    #     return protein_infos

    def get_dna(self):
        try:
            if not self.coding_sequence:
                return

            return DNARecord.get_id(self.coding_sequence[0].id)

        except Exception as e:
            print("The DNA sequence could not be retrieved. The error is: ", e)
            return

    def _nblast(sequence: str, n_hits: int = None) -> List["ProteinRecord"]:
        # blast_record = NCBIXML.read(result_handle)
        raise NotImplementedError("This method is not implemented yet.")

    def fasta_string(self) -> str:
        return f">{self.id}\n{self.sequence}"

    def to_fasta(self, path: str):
        with open(path, "w") as f:
            f.write(self.fasta_string())

        print(f"💾 Sequence saved to {path}")

    @staticmethod
    def _get_accessions(blast_record: NCBIXML) -> List[str]:
        return [alignment.accession for alignment in blast_record.alignments]

    def from_ncbi(self):
        raise DeprecationWarning("This method is deprecated. Use `get_id` instead.")

    def from_accessions(self):
        raise DeprecationWarning("This method is deprecated. Use `get_ids` instead.")
