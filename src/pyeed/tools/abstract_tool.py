import shutil
import tempfile
from abc import ABC, abstractmethod
from enum import Enum
from typing import Any

import httpx
from pydantic import BaseModel, PrivateAttr


class ServiceURL(Enum):
    CLUSTALO = "http://clustalo:5001/align"
    BLAST = "http://blast:6001/"
    FOLDSEEK = "http://foldseek:7001/foldseek"
    MMSEQS = "http://mmseqs:8001/easycluster"


class AbstractTool(BaseModel, ABC):
    """
    Abstract base class for tools.
    """

    _service_url: str = PrivateAttr()
    _tempdir_path: str = PrivateAttr()

    def model_post_init(self, __context: Any) -> None:
        """Create temporary directory."""
        self._tempdir_path = tempfile.mkdtemp()

    def _delete_temp_dir(self) -> None:
        """Deletes the temporary directory."""
        shutil.rmtree(self._tempdir_path)

    @abstractmethod
    def run_service(self, data: dict[str, Any]) -> httpx.Response:
        """Executes the service."""
        pass
