import logging
import secrets
import logging.config
from pathlib import Path
from abc import ABC, abstractmethod
from typing import Any


# path_config = Path(__file__).parent.parent.parent / "logging.conf"
# logging.config.fileConfig(path_config)
# LOGGER = logging.getLogger("pyeed")

LOGGER = logging.getLogger(__name__)


class AbstractFetcher(ABC):
    def __init__(self, foreign_id: str):
        super().__init__()
        self.foreign_id = foreign_id

    @abstractmethod
    def get(self):
        pass

    @abstractmethod
    def map(self, handle: Any, cls):
        pass

    @abstractmethod
    def fetch(self):
        pass

    @staticmethod
    def get_substitute_email() -> str:
        return f"{secrets.token_hex(8)}@gmail.com"
