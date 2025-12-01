from .utils.environment import IN_IPYTHON, IN_NOTEBOOK
from .utils.logging_setup import setup_logging

setup_logging()

__all__ = [
    "IN_IPYTHON",
    "IN_NOTEBOOK",
]
