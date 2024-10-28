import logging
import logging.config
import shutil
import httpx

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
    def run_service(self, data) -> httpx.Response:
        """Executes the service."""
        pass
