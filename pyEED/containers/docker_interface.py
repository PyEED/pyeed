import os
import shutil
from typing import Any, List
import tempfile
from abc import ABC, abstractmethod
from pydantic import BaseModel

import docker
from docker.models.containers import Container
from docker.models.images import Image
from docker.client import DockerClient


from pyEED.containers.containers import ToolContainer
from pyEED.core.abstractsequence import AbstractSequence
from pyEED.utility.utilities import create_multifaster


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
        get_image(): Gets the image from Docker Hub.
        run_container(command: str, data: Any): Runs a container and executes a command.
        setup_input_data(data: Any) -> str: Creates the input data for the container.
        extract_output_data(): Extracts the output data from the container.
        setup_command() -> str: Creates the command to be executed in the container.
        delete_temp_dir(): Deletes the temporary directory.
    """

    _container_info: ToolContainer
    _client: DockerClient = docker.from_env()
    _tempdir_path: str = tempfile.mkdtemp()

    def get_image(self) -> Image:
        """Gets the image from Docker Hub. If the image is not found, it will be pulled."""
        try:
            return self._client.images.get(self._container_info.value)
        except docker.errors.ImageNotFound:
            return self._client.images.pull(self._container_info.value)

    def run_container(self, command: str, data: Any) -> Container:
        """Runs a container from a given image and executes a command.
        Data is mounted to the container."""
        image = self.get_image()
        self.setup_input_data(data=data)

        container = self._client.containers.run(
            image=image,
            command=command,
            name=self._container_info.name,
            auto_remove=True,
            volumes={self._tempdir_path: {"bind": "/data/", "mode": "rw"}},
        )

        return container

    def delete_temp_dir(self):
        """Deletes the temporary directory."""
        shutil.rmtree(self._tempdir_path)

    @abstractmethod
    def setup_input_data(self, data: Any) -> str:
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


class ClustalOmega(AbstractContainer):
    _container_info: ToolContainer = ToolContainer.CLUSTALO

    def setup_input_data(self, data: List[AbstractSequence]) -> str:
        fasta = create_multifaster(data)
        path = os.path.join(self._tempdir_path, "input.fasta")

        with open(path, "w") as file:
            file.write(fasta)

    def setup_command(self):
        return "clustalo -i /data/input.fasta -o /data/output.fasta"

    def extract_output_data(self):
        import Bio.AlignIO

        with open(os.path.join(self._tempdir_path, "output.fasta"), "r") as file:
            alignment = Bio.AlignIO.read(file, "fasta")

        self.delete_temp_dir()
        return alignment

    def align(self, sequences: List[AbstractSequence]):
        """Aligns multiple sequences and returns the alignment result.

        Args:
            sequences (List[AbstractSequence]): List of sequences to be aligned.

        Returns:
            _type_: _description_
        """
        self.run_container(
            command=self.setup_command(),
            data=sequences,
        )
        return self.extract_output_data()
