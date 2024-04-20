import asyncio
import os
import warnings
from concurrent.futures import ThreadPoolExecutor
from typing import List, Optional

from IPython.display import clear_output
from pydantic import Field
from rich.console import Console
from rich.status import Status
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature

from pyeed.container.abstract_container import Blastp
from pyeed.core.site import Site


@forge_signature
class ProteinInfo(AbstractSequence):
    """Description of a protein sequence. Additionally, the `ProteinSequence` contains annotations for sites and regions of the protein sequence alongside information on the organism. Furthermore, the `ProteinSequence` contains information on the coding sequence of the protein sequence, which allows later retrieval of the corresponding nucleotide sequence."""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("proteininfoINDEX"),
        xml="@id",
    )

    family_name: Optional[str] = Field(
        default=None,
        description="Family name of the protein",
    )

    regions: List[ProteinRegion] = Field(
        description="Domains of the protein",
        default_factory=ListPlus,
        multiple=True,
    )

    sites: List[Site] = Field(
        description="Annotations of different sites",
        default_factory=ListPlus,
        multiple=True,
    )

    coding_sequence_ref: Optional[DNARegion] = Field(
        description="Defines the coding sequence of the protein",
        default_factory=DNARegion,
    )

    ec_number: Optional[str] = Field(
        default=None,
        description="Enzyme Commission number",
    )

    mol_weight: Optional[float] = Field(
        default=None,
        description="Calculated molecular weight of the protein",
    )

    substrates: List[Substrate] = Field(
        description="Promiscuous substrates of the protein",
        default_factory=ListPlus,
        multiple=True,
    )

    def add_to_regions(
        self,
        type: Optional[ProteinRegionType] = None,
        name: Optional[str] = None,
        spans: List[Span] = ListPlus(),
        note: Optional[str] = None,
        cross_reference: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'ProteinRegion' to attribute regions

        Args:
            id (str): Unique identifier of the 'ProteinRegion' object. Defaults to 'None'.
            type (): Type of the region within the protein sequence. Defaults to None
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
        self.regions.append(ProteinRegion(**params))
        return self.regions[-1]

    def add_to_sites(
        self,
        name: Optional[str] = None,
        type: Optional[ProteinSiteType] = None,
        positions: List[int] = ListPlus(),
        cross_ref: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Site' to attribute sites

        Args:
            id (str): Unique identifier of the 'Site' object. Defaults to 'None'.
            name (): Name of the site. Defaults to None
            type (): Type of the site. Defaults to None
            positions (): Positions of the site. Defaults to ListPlus()
            cross_ref (): Database cross reference. Defaults to None
        """
        params = {
            "name": name,
            "type": type,
            "positions": positions,
            "cross_ref": cross_ref,
        }
        if id is not None:
            params["id"] = id
        self.sites.append(Site(**params))
        return self.sites[-1]

    def add_to_substrates(
        self,
        name: Optional[str] = None,
        inchi: Optional[str] = None,
        smiles: Optional[str] = None,
        chebi_id: Optional[str] = None,
        citation: Optional[Citation] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Substrate' to attribute substrates

        Args:
            id (str): Unique identifier of the 'Substrate' object. Defaults to 'None'.
            name (): Name of the substrate. Defaults to None
            inchi (): InChI code of the substrate. Defaults to None
            smiles (): SMILES code of the substrate. Defaults to None
            chebi_id (): ChEBI ID of the substrate. Defaults to None
            citation (): Citations of the substrate. Defaults to None
        """
        params = {
            "name": name,
            "inchi": inchi,
            "smiles": smiles,
            "chebi_id": chebi_id,
            "citation": citation,
        }
        if id is not None:
            params["id"] = id
        self.substrates.append(Substrate(**params))
        return self.substrates[-1]

    @classmethod
    def get_id(cls, protein_id: str) -> "ProteinInfo":
        import nest_asyncio

        from pyeed.fetch.proteinfetcher import ProteinFetcher

        nest_asyncio.apply()

        """
        This method creates a 'ProteinInfo' object from a given protein accession ID.

        Args:
            protein_id (str): ID of the protein in NCBI or UniProt database.

        Returns:
            ProteinInfo: 'ProteinInfo' with information of the corresponding protein_id.
        """

        if isinstance(protein_id, list) and all(isinstance(x, str) for x in protein_id):
            warnings.warn("For getting multiple sequences by ID use `get_ids` instead.")
            return cls.get_ids(protein_id)

        sequences = asyncio.run(ProteinFetcher(ids=[protein_id]).fetch(quiet=True))[0]
        clear_output()
        return sequences

    @classmethod
    def get_ids(cls, accession_ids: List[str]) -> List["ProteinInfo"]:
        import nest_asyncio

        from pyeed.fetch.proteinfetcher import ProteinFetcher

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

        import nest_asyncio

        from pyeed.fetch.blast import Blast, BlastProgram, NCBIDataBase
        from pyeed.fetch.proteinfetcher import ProteinFetcher

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
        database: str = "nr",
        matrix: str = "BLOSUM62",
        identity: float = 0.0,
        **kwargs,
    ) -> List["ProteinInfo"]:
        """
        Runs a BLAST search using the NCBI BLAST service to find similar protein sequences.

        Args:
            n_hits (int): The number of hits to retrieve.
            e_value (float, optional): The maximum E-value threshold for reporting hits. Defaults to 10.0.
            database (str, optional): The database to search against. Defaults to "nr".
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

        from pyeed.fetch.blast import Blast, BlastProgram, NCBIDataBase
        from pyeed.fetch.proteinfetcher import ProteinFetcher

        nest_asyncio.apply()

        assert database in NCBIDataBase

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
            result = asyncio.run(blaster.async_run(database, program, executor))
            clear_output()

        accessions = blaster.extract_accession(result)

        return asyncio.run(ProteinFetcher(ids=accessions).fetch(force_terminal=False))

    def blastp(
        self,
        db_path: str,
        identity: float = 0,
        evalue: float = 10,
        n_hits: int = 500,
        subst_matrix: str = "BLOSUM62",
        word_size: int = 3,
        gapopen: int = 11,
        gapextend: int = 1,
        threshold: int = 11,
        n_cores: int = os.cpu_count(),
        ncbi_key: str = None,
        email: str = None,
    ):
        blaster = Blastp(
            _db_path=db_path,
            identity=identity,
            evalue=evalue,
            n_hits=n_hits,
            subst_matrix=subst_matrix,
            word_size=word_size,
            gapopen=gapopen,
            gapextend=gapextend,
            threshold=threshold,
            n_cores=n_cores,
            ncbi_key=ncbi_key,
        )

        command = blaster.setup_command()
        accession_ids = blaster.run_container(
            command=command, data=self._fasta_string()
        )
        protein_infos = ProteinInfo.get_ids(accession_ids, email, ncbi_key)
        protein_infos.insert(0, self)
        return protein_infos

    def get_dna(self):
        if not self.coding_sequence_ref:
            return

        return DNAInfo.from_ncbi(self.coding_sequence_ref.id)

    def _nblast(sequence: str, n_hits: int = None) -> List["ProteinInfo"]:
        # blast_record = NCBIXML.read(result_handle)
        raise NotImplementedError("This method is not implemented yet.")

    @staticmethod
    def _get_accessions(blast_record: NCBIXML) -> List[str]:
        return [alignment.accession for alignment in blast_record.alignments]

    def from_ncbi(self):
        raise DeprecationWarning("This method is deprecated. Use `get_id` instead.")

    def from_accessions(self):
        raise DeprecationWarning("This method is deprecated. Use `get_ids` instead.")


if __name__ == "__main__":
    seq_string = "MSDRNIRVEPVVGRAVEEQDVEIVERKGLGHPDSLCDGIAEHVSQALARAYIDRVGKVLHYNTDETQLVAGTAAPAFGGGEVVDPIYLLITGRATKEYEGTKIPAETIALRAAREYINETLPFLEFGTDVVVDVKLGEGSGDLQEVFGEDGKQVPMSNDTSFGVGHAPLTETERIVLEAERALNGDYSDDNPAVGQDIKVMGKREGDDIDVTVAVAMVDRYVDDLDGYEAAVAGVREFVADLATDYTDRNVSVHVNTADDYDEGAIYLTTTGTSAEQGDDGSVGRGNRSNGLITPNRSMSMEATSGKNPVNHIGKIYNLLSTEIARTVVDEVDGIREIRIRLLSQIGQPIDKPHVADANLVTEDGIEIADIEDEVEAIIDAELENVTSITERVIDGELTTF"

    seq = ProteinInfo.from_sequence(seq_string)
    print(seq)
