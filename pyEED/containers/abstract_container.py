import docker
import shutil
from enum import Enum
from typing import Any
import tempfile
from abc import ABC, abstractmethod
from pydantic import BaseModel, PrivateAttr

from docker.client import DockerClient
from docker.models.containers import Container
from docker.models.images import Image


class ToolImage(Enum):
    CLUSTALO = "biocontainers/clustal-omega:v1.2.1_cv5"


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
        self._client = docker.from_env()

    def get_image(self) -> Image:
        """Gets the image from Docker Hub. If the image is not found, it will be pulled."""
        try:
            return self._client.images.get(self._container_info.value)
        except docker.errors.ImageNotFound:
            print(f"⬇️ Pulling {self._container_info.name}")
            return self._client.images.pull(self._container_info.value)

    def run_container(self, command: str, data: Any) -> Container:
        """Runs a container from a given image and executes a command.
        Data is mounted to the container."""
        try:
            image = self.get_image()
            self.create_file(data=data)

            print(f"🏃 Running {self._container_info.name}")
            self._client.containers.run(
                image=image,
                command=command,
                name=self._container_info.name,
                auto_remove=True,
                volumes={self._tempdir_path: {"bind": "/data/", "mode": "rw"}},
            )
        except Exception as e:
            print(f"Error running container: {e}")

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
