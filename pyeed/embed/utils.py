import logging
import os
from contextlib import contextmanager

import torch
from dotenv import load_dotenv
from huggingface_hub import login
from loguru import logger


def _login_hf() -> str:
    """Login to Hugging Face."""
    load_dotenv()
    token = os.getenv("HUGGINGFACE_HUB_TOKEN")
    if isinstance(token, str):
        login(token=token)
        return token
    else:
        raise ValueError("HUGGINGFACE_HUB_TOKEN environment variable is not set")


def _free_device_memory(device_ids: list[int] | None = None) -> None:
    """
    Free GPU memory on specified devices.

    Args:
        device_ids: List of GPU device IDs to free memory on.
                   If None, frees memory on all available devices.
    """
    if not torch.cuda.is_available():
        return

    if device_ids is None:
        device_ids = list(range(torch.cuda.device_count()))

    for device_id in device_ids:
        with torch.cuda.device(device_id):
            torch.cuda.empty_cache()
            torch.cuda.synchronize()

    logger.debug(
        "Freed GPU memory on devices",
        extra={
            "device_ids": device_ids,
        },
    )


_SILENCE_SUBSTRINGS = (
    "were not initialized from the model checkpoint",
    "You should probably TRAIN this model on a down-stream task",
)


class _TransformersInitFilter(logging.Filter):
    """Filter to silence Transformers initialization warnings."""

    def filter(self, record: logging.LogRecord) -> bool:
        msg = record.getMessage()
        return not any(s in msg for s in _SILENCE_SUBSTRINGS)


@contextmanager
def silence_transformers_init_only():
    """Context manager to silence Transformers initialization warnings."""
    logger = logging.getLogger("transformers.modeling_utils")
    flt = _TransformersInitFilter()
    logger.addFilter(flt)
    # ensure WARNINGs flow (we only drop the exact ones)
    prev_level = logger.level
    logger.setLevel(min(prev_level or logging.WARNING, logging.WARNING))
    try:
        yield
    finally:
        logger.removeFilter(flt)
        logger.setLevel(prev_level)
