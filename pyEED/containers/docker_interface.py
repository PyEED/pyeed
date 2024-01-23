from typing import List
import tempfile
from abc import ABC, abstractmethod
from pydantic import BaseModel

import docker
from docker.models.containers import Container
from docker.models.images import Image


from pyEED.containers.containers import ToolContainer
from pyEED.core.abstractsequence import AbstractSequence
from pyEED.utility.utilities import create_multifaster


class AbstractContainer(BaseModel, ABC):
    _container_info: ToolContainer
    _client = docker.from_env()

    def get_image(self) -> Image:
        """Gets the image from dockerhub. If the image is not found, it will be pulled."""
        try:
            return self._client.images.get(self._container_info.value)
        except docker.errors.ImageNotFound:
            return self._client.images.pull(self._container_info.value)

    def run_container(self, command: str) -> Container:
        """Gets the container and runs it."""
        image = self.get_image()
        container = self._client.containers.run(
            image=image,
            command=command,
            name=self._container_info.name,
            auto_remove=True,
            volumes={self.mount_data().name: {"bind": "temp/", "mode": "rw"}},
        )

        return container

    @abstractmethod
    def construct_command(self):
        pass

    @abstractmethod
    def mount_data(self):
        pass


class ClustalOmega(AbstractContainer):
    _container_info: ToolContainer = ToolContainer.CLUSTALO

    def mount_data(self, sequences: List[AbstractSequence]):
        print(create_multifaster(sequences))
        with tempfile.NamedTemporaryFile(mode="w+") as temp:
            temp.write(create_multifaster(sequences))
            temp.flush()
            print(temp.name)

            return temp

    def construct_command(self):
        pass

    def align(self, sequences: List[AbstractSequence]):
        command = self.construct_command()


if __name__ == "__main__":
    clustal = ClustalOmega()
    c = clustal.run_container("clustalo --help")

    print(c)
