from __future__ import annotations

import asyncio
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Protocol

from rich.progress import Progress, TaskID

SENTINEL = object()


@dataclass
class PipelineContext:
    """Shared state across pipeline stages."""

    progress: Progress
    added_nodes: defaultdict[str, set[str]] = field(default_factory=lambda: defaultdict(set))
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
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        """Run the stage until SENTINEL received on all input queues.

        Args:
            input_queues: Named input queues (empty dict for source stages)
            output_queues: Named output queues (empty dict for sink stages)
            context: Shared pipeline context
            progress: Progress object for reporting progress
            task_id: TaskID for progress tracking
        """
        ...
