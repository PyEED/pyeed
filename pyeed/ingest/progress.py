from __future__ import annotations

from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from ..environment import IN_NOTEBOOK

__all__ = ["create_progress"]

CONSOLE = Console(force_jupyter=IN_NOTEBOOK, force_terminal=not IN_NOTEBOOK)


def create_progress(transient: bool = False, progress: Progress | None = None) -> Progress:
    """Create or reuse a Rich progress bar with standard configuration.

    Args:
        transient: If True, progress bar disappears after completion.
        progress: If provided, returns the existing Progress instance instead of creating a new one.
                  This allows sharing a single progress context across multiple operations.

    Returns:
        Configured Progress instance (either new or reused).
    """
    if progress is not None:
        return progress

    return Progress(
        SpinnerColumn(),
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        transient=transient,
        console=CONSOLE,
    )
