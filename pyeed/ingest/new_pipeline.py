from __future__ import annotations

import asyncio
from dataclasses import dataclass, field
from typing import Protocol

from rich.progress import Progress, TaskID

from pyeed.utils.progress import ProgressReporter

from ..utils.progress import create_progress


class PipelineStage(Protocol):
    """Protocol for pipeline stages. Stages should implement run() method."""

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue],
        output_queues: dict[str, asyncio.Queue],
        progress_reporters: dict[str, ProgressReporter],
        context: PipelineContext,
    ) -> None:
        """Run the stage until SENTINEL received on all inputs."""
        ...


@dataclass
class PipelineContext:
    """Shared state across pipeline stages."""

    progress: Progress
    skip_neo4j: set[str] = field(default_factory=set)
    skip_milvus: set[str] = field(default_factory=set)
    stats: dict[str, int] = field(default_factory=dict)


@dataclass
class StageConfig:
    """Configuration for a pipeline stage."""

    stage: PipelineStage
    input_queues: list[str]
    output_queues: list[str]
    progress_tasks: dict[str, str]


class Pipeline:
    def __init__(
        self,
        progress: Progress | None = None,
    ):
        self.stages: list[StageConfig] = []
        self.queues: dict[str, asyncio.Queue] = {}
        self.progress = progress or create_progress()
        self.context = PipelineContext(progress=self.progress)
        self._progress_tasks: dict[str, TaskID] = {}

    def add_queue(self, name: str, maxsize: int) -> str:
        """Register a named queue."""
        self.queues[name] = asyncio.Queue(maxsize=maxsize)
        return name

    def add_stage(
        self,
        stage: PipelineStage,
        input_queues: list[str],
        output_queues: list[str],
        progress_tasks: dict[str, str] | None = None,
    ) -> None:
        """Register a stage with queues and progress tasks.

        Args:
            stage: Stage implementation
            input_queues: List of input queue names
            output_queues: List of output queue names
            progress_tasks: Dict mapping progress names to descriptions
                           e.g., {"read": "Read FASTA", "embed": "Embed"}
        """
        config = StageConfig(
            stage=stage,
            input_queues=input_queues,
            output_queues=output_queues,
            progress_tasks=progress_tasks or {},
        )
        self.stages.append(config)

    async def run(
        self,
        total_items: dict[str, int] | None = None,
    ) -> dict[str, int]:
        """Run all stages concurrently.

        Args:
            total_items: Dict mapping progress task names to totals
                        e.g., {"read": 10000, "embed": 10000}
        """
        # Create all progress tasks upfront
        for config in self.stages:
            for task_name, task_desc in config.progress_tasks.items():
                if task_name not in self._progress_tasks:
                    total = (total_items or {}).get(task_name)
                    task_id = self.progress.add_task(task_desc, total=total)
                    self._progress_tasks[task_name] = task_id

        # Launch all stages
        tasks = []
        for config in self.stages:
            in_qs = {name: self.queues[name] for name in config.input_queues}
            out_qs = {name: self.queues[name] for name in config.output_queues}

            # Build ProgressReporters for this stage
            reporters = {
                name: ProgressReporter(self.progress, self._progress_tasks[name])
                for name in config.progress_tasks.keys()
            }

            tasks.append(
                asyncio.create_task(config.stage.run(in_qs, out_qs, reporters, self.context))
            )

        await asyncio.gather(*tasks)
        return self.context.stats
