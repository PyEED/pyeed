from __future__ import annotations

import asyncio
from dataclasses import dataclass, field
from typing import Any, Protocol

from rich.progress import Progress

from ...utils.progress import ProgressReporter

SENTINEL = object()


@dataclass
class PipelineContext:
    """Shared context across all pipeline stages."""

    progress: Progress
    skip_neo4j: set[str] = field(default_factory=set)
    skip_milvus: set[str] = field(default_factory=set)
    stats: dict[str, int] = field(default_factory=dict)


class PipelineStage(Protocol):
    """Protocol for pipeline stages.

    Stages consume from input queue(s), process data, and produce to output queue(s).
    Progress reporting is optional but recommended for observability.
    """

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        progress_reporters: dict[str, ProgressReporter],
        context: PipelineContext,
    ) -> None:
        """Run the stage until SENTINEL received on all input queues.

        Args:
            input_queues: Named input queues (empty dict for source stages)
            output_queues: Named output queues (empty dict for sink stages)
            progress_reporters: Named progress reporters for this stage
            context: Shared pipeline context
        """
        ...
