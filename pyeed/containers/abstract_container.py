import os
import logging
import logging.config
import docker
import shutil
import random
import string
from pathlib import Path
from enum import Enum
from typing import Any
import tempfile
from abc import ABC, abstractmethod
from pydantic import BaseModel, Field, PrivateAttr

from docker.client import DockerClient
from docker.models.containers import Container
from docker.models.images import Image
from docker.errors import DockerException

path_config = Path(__file__).parent.parent.parent / "logging.conf"
logging.config.fileConfig(path_config)
logger = logging.getLogger("pyeed")


class ToolImage(Enum):
    CLUSTALO = "biocontainers/clustal-omega:v1.2.1_cv5"
    BLAST = "haeussma/blast:latest"
    MMSEQS2 = "soedinglab/mmseqs2:master"


class AbstractContainer(BaseModel, ABC):
    """
    Abstract base class for containers.

    This class provides a common interface for interacting with Docker containers.
    Subclasses should implement the `setup_input_data`, `extract_output_data` and
    `setup_command` methods.

    Attributes:
        _container_info (ToolContainer): Information about the container.
        _client (DockerClient): Docker client instance.
        _tempdir_path (str): Path to the temporary directory.

    Methods:
        get_image() -> Image:
            Gets the image from Docker Hub. If the image is not found, it will be pulled.

        run_container(command: str, data: Any) -> Container:
            Runs a container from a given image and executes a command.
            Data is mounted to the container.

        setup_input_data(data: Any) -> str:
            Creates the input data for the container.

        extract_output_data():
            Extracts the output data from the container.

        setup_command() -> str:
            Creates the command to be executed in the container.

        _delete_temp_dir():
            Deletes the temporary directory.
    """

    _container_info: ToolImage

    _client: DockerClient = PrivateAttr()
    _tempdir_path: str = PrivateAttr()

    def __init__(self, **kwargs):
        BaseModel.__init__(self, **kwargs)
        super().__init__(**kwargs)
        self._tempdir_path = tempfile.mkdtemp()
        self._client = self._initialize_docker_client()
        logger.debug(
            f"Successfully initialized container of type: {self._container_info}, _tempdir_path: {self._tempdir_path}"
        )

    def _initialize_docker_client(self) -> DockerClient:

        logger.debug("Initializing Docker client")
        try:
            client = docker.from_env()
            client.ping()
            return client
        except DockerException as e:
            logger.error(f"Docker is not running. Start the Docker application. {e}")

    def get_image(self) -> Image:
        """Gets the image from Docker Hub. If the image is not found, it will be pulled."""
        try:
            logger.debug(f"Getting image {self._container_info.value}")
            return self._client.images.get(self._container_info.value)
        except docker.errors.ImageNotFound:
            logger.info(f"⬇️ Pulling {self._container_info.name}")

            return self._client.images.pull(self._container_info.value)

    def run_container(self, command: str, data: Any) -> Container:
        """Runs a container from a given image and executes a command.
        Data is mounted to the container."""
        try:
            image = self.get_image()
            self.create_file(data=data)

            logger.info(f"🏃 Running {self._container_info.name}")
            self._client.containers.run(
                image=image,
                command=command,
                name=self._container_info.name
                + "_"
                + "".join(random.choices(string.ascii_lowercase, k=5)),
                auto_remove=True,
                volumes={self._tempdir_path: {"bind": "/data/", "mode": "rw"}},
            )
        except Exception as e:
            logger.error(f"Error running {self._container_info} container: {e}")

    def _delete_temp_dir(self):
        """Deletes the temporary directory."""
        shutil.rmtree(self._tempdir_path)

    @abstractmethod
    def create_file(self, data: Any) -> str:
        """Creates the input data for the container."""
        pass

    @abstractmethod
    def extract_output_data(self):
        """Extracts the output data from the container."""
        pass

    @abstractmethod
    def setup_command(self) -> str:
        """Creates the command which is executed in the container."""
        pass


class Blastp(AbstractContainer):
    """
    Class for running BLASTP.

    Attributes:
        _container_info (ToolContainer): Information about the container.
        _client (DockerClient): Docker client instance.
        _tempdir_path (str): Path to the temporary directory.

    Methods:
        create_file(data: Any) -> str:
            Creates the input data for the container.

        extract_output_data():
            Extracts the output data from the container.

        setup_command() -> str:
            Creates the command to be executed in the container.
    """

    identity: float = Field(default=0.0, description="Minimum identity to safe hits.")
    evalue: float = Field(default=10, description="Expectation value (E) to safe hits.")
    n_hits: int = Field(default=500, description="Maximum number of hits to return.")
    subs_matrix: str = Field(default="BLOSUM62")
    word_size: int = Field(
        default=3, ge=2, le=7, description="Word size of the initial match."
    )
    gapopen: int = Field(default=11, description="Cost to open a gap.")
    gapextend: int = Field(default=1, description="Cost to extend a gap.")
    threshold: int = Field(
        default=11, description="Minimum score to add a word to the BLAST lookup table."
    )

    _container_info = ToolImage.BLAST
    _db_path: str = PrivateAttr()
    _n_cores: int = PrivateAttr(default=os.cpu_count())
    _ncbi_key: str = PrivateAttr(default=None)

    def __init__(self, _db_path: str, ncbi_key: str = None, **kwargs):
        super().__init__(**kwargs)
        self._db_path = _db_path
        if not self._ncbi_key:
            try:
                self._ncbi_key = os.environ["NCBI_API_KEY"]
            except KeyError:
                self._ncbi_key = ncbi_key

    def run_container(self, command: str, data: Any) -> Container:
        try:
            image = self.get_image()
            self.create_file(data=data)

            print(f"🏃 Running {self._container_info.name}")
            self._client.containers.run(
                image=image,
                command=command,
                name=self._container_info.name
                + "".join(random.choices(string.ascii_lowercase, k=5)),
                auto_remove=True,
                volumes={
                    self._tempdir_path: {"bind": "/data/", "mode": "rw"},
                    self._db_path: {"bind": "/db/", "mode": "rw"},
                },
            )
            return self.extract_output_data()
        except Exception as e:
            print(f"Error running container: {e}")

    def create_file(self, data: str) -> str:
        """Creates the input data for the container."""

        with open(f"{self._tempdir_path}/input.fasta", "w") as f:
            f.write(data)
            logger.debug(f"Created input file: {self._tempdir_path}/blastp.fasta")
        return f"{self._tempdir_path}/blastp.fasta"

    def extract_output_data(self):
        from Bio import SearchIO

        """Extracts the output data from the container."""

        search_io = SearchIO.read(f"{self._tempdir_path}/blastp.out", "blast-tab")
        logger.debug("Extracted output data")

        # TODO add debug info on how many sequences were discarded due to identity and evalue thresholds

        return [
            result.id
            for result in search_io
            if result._items[0].ident_pct >= self.identity * 100
        ]

    def setup_command(self) -> str:
        """Creates the command which is executed in the container."""
        logger.debug(f"Setting up command for {self._container_info}")
        return (
            f"blastp "
            f"-query /data/input.fasta "
            f"-db /db/clustered_db "
            f"-out /data/blastp.out "
            f"-outfmt 6 "
            f"-gapopen {self.gapopen} "
            f"-gapextend {self.gapextend} "
            f"-threshold {self.threshold} "
            f"-num_threads {self._n_cores} "
            f"-matrix {self.subs_matrix} "
            f"-word_size {self.word_size} "
            f"-evalue {self.evalue} "
            f"-max_target_seqs {self.n_hits}"
        )
