import logging
import os

LOG_FILE = os.getenv("PYEED_LOG_FILE", "pyeed.log")
LOG_LEVEL = os.getenv("PYEED_LOG_LEVEL", "INFO").upper()


def setup_logging() -> None:
    """Configure global logging for the whole package."""
    logging.basicConfig(
        level=LOG_LEVEL,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        filename=LOG_FILE,
        filemode="a",  # append
    )
