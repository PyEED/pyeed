import shutil
import tempfile
from abc import ABC, abstractmethod
from enum import Enum

import httpx
from pydantic import BaseModel, PrivateAttr


class ServiceURL(Enum):
    CLUSTALO = "http://clustalo:5001/align"
    BLASTP = "http://blast:6001/blastp"
    BLASTN = "http://blast:6001/blastn"
    FOLDSEEK = "http://foldseek:7001/foldseek"
    MMSEQS = "http://mmseqs:8001/easycluster"


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
    def run_service(self, data) -> httpx.Response:
        """Executes the service."""
        pass
