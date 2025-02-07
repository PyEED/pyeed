import os
from typing import Generator

import pytest


@pytest.fixture(autouse=True)
def skip_hf_login() -> Generator[None, None, None]:
    """Skip Hugging Face login during tests."""
    os.environ["PYTEST_DISABLE_HF_LOGIN"] = "1"
    yield
    if "PYTEST_DISABLE_HF_LOGIN" in os.environ:
        del os.environ["PYTEST_DISABLE_HF_LOGIN"]
