# yourpkg/env.py
from __future__ import annotations


def in_notebook() -> bool:
    """
    True if running in a Jupyter-like kernel (classic, JLab, VSCode, Colab).
    Safe if IPython isn't installed.
    """
    try:
        # If an IPython shell is present, check its type
        from IPython import get_ipython  # type: ignore

        ip = get_ipython()
        if ip is None:
            return False
        # ZMQInteractiveShell => Jupyter/VSCode/Colab kernels
        return ip.__class__.__name__ == "ZMQInteractiveShell"
    except Exception:
        return False


def in_ipython() -> bool:
    try:
        from IPython import get_ipython  # type: ignore

        return get_ipython() is not None
    except Exception:
        return False


# Cached flags (computed once at import)
IN_NOTEBOOK: bool = in_notebook()
IN_IPYTHON: bool = in_ipython()
