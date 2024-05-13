import logging
import logging.config
import shutil
import tempfile
from abc import ABC, abstractmethod
from enum import Enum

from pydantic import BaseModel, PrivateAttr

# path_config = Path(__file__).parent.parent.parent / "logging.conf"
# logging.config.fileConfig(path_config)
# logger = logging.getLogger("pyeed")

logger = logging.getLogger(__name__)


class ServiceURL(Enum):
    CLUSTALO = "http://clustalo:5001/align"


class AbstractTool(BaseModel, ABC):
    """
    Abstract base class for tools.
    """

    _service_url: str = PrivateAttr()
    _tempdir_path: str = PrivateAttr()

    def model_post_init(self, __context) -> None:
        """Create temporary directory."""
        self._tempdir_path = tempfile.mkdtemp()

    def _delete_temp_dir(self):
        """Deletes the temporary directory."""
        shutil.rmtree(self._tempdir_path)

    @abstractmethod
    def run_service(self):
        """Executes the service."""
        pass


# class Blastp(AbstractContainer):
#     """
#     Class for running BLASTP.

#     Attributes:
#         _container_info (ToolContainer): Information about the container.
#         _client (DockerClient): Docker client instance.
#         _tempdir_path (str): Path to the temporary directory.

#     Methods:
#         create_file(data: Any) -> str:
#             Creates the input data for the container.

#         extract_output_data():
#             Extracts the output data from the container.

#         setup_command() -> str:
#             Creates the command to be executed in the container.
#     """

#     identity: float = Field(default=0.0, description="Minimum identity to safe hits.")
#     evalue: float = Field(default=10, description="Expectation value (E) to safe hits.")
#     n_hits: int = Field(default=500, description="Maximum number of hits to return.")
#     subs_matrix: str = Field(default="BLOSUM62")
#     word_size: int = Field(
#         default=3, ge=2, le=7, description="Word size of the initial match."
#     )
#     gapopen: int = Field(default=11, description="Cost to open a gap.")
#     gapextend: int = Field(default=1, description="Cost to extend a gap.")
#     threshold: int = Field(
#         default=11, description="Minimum score to add a word to the BLAST lookup table."
#     )

#     _container_info = ToolImage.BLAST
#     _db_path: str = PrivateAttr()
#     _n_cores: int = PrivateAttr(default=os.cpu_count())
#     _ncbi_key: str = PrivateAttr(default=None)

#     def __init__(self, _db_path: str, ncbi_key: str = None, **kwargs):
#         super().__init__(**kwargs)
#         self._db_path = _db_path
#         if not self._ncbi_key:
#             try:
#                 self._ncbi_key = os.environ["NCBI_API_KEY"]
#             except KeyError:
#                 self._ncbi_key = ncbi_key

#     def run_container(self, command: str, data: Any) -> Container:
#         try:
#             image = self.get_image()
#             self.create_file(data=data)

#             print(f"🏃 Running {self._container_info.name}")
#             self._client.containers.run(
#                 image=image,
#                 command=command,
#                 name=self._container_info.name
#                 + "".join(random.choices(string.ascii_lowercase, k=5)),
#                 auto_remove=True,
#                 volumes={
#                     self._tempdir_path: {"bind": "/data/", "mode": "rw"},
#                     self._db_path: {"bind": "/db/", "mode": "rw"},
#                 },
#             )
#             return self.extract_output_data()
#         except Exception as e:
#             print(f"Error running container: {e}")

#     def create_file(self, data: str) -> str:
#         """Creates the input data for the container."""

#         with open(f"{self._tempdir_path}/input.fasta", "w") as f:
#             f.write(data)
#             logger.debug(f"Created input file: {self._tempdir_path}/blastp.fasta")
#         return f"{self._tempdir_path}/blastp.fasta"

#     def extract_output_data(self):
#         from Bio import SearchIO

#         """Extracts the output data from the container."""

#         search_io = SearchIO.read(f"{self._tempdir_path}/blastp.out", "blast-tab")
#         logger.debug("Extracted output data")

#         # TODO add debug info on how many sequences were discarded due to identity and evalue thresholds

#         return [
#             result.id
#             for result in search_io
#             if result._items[0].ident_pct >= self.identity * 100
#         ]

#     def setup_command(self) -> str:
#         """Creates the command which is executed in the container."""
#         logger.debug(f"Setting up command for {self._container_info}")
#         return (
#             f"blastp "
#             f"-query /data/input.fasta "
#             f"-db /db/clustered_db "
#             f"-out /data/blastp.out "
#             f"-outfmt 6 "
#             f"-gapopen {self.gapopen} "
#             f"-gapextend {self.gapextend} "
#             f"-threshold {self.threshold} "
#             f"-num_threads {self._n_cores} "
#             f"-matrix {self.subs_matrix} "
#             f"-word_size {self.word_size} "
#             f"-evalue {self.evalue} "
#             f"-max_target_seqs {self.n_hits}"
#         )
