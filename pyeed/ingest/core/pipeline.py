from __future__ import annotations

import asyncio
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from rich.progress import Progress, TaskID

from ...utils.progress import create_progress
from ..model.pyeedbase import PyeedBase
from .protocol import PipelineContext, PipelineStage


@dataclass
class ChildRecord[T: PyeedBase]:
    data: list[T]
    parent_label: str
    parent_field: str
    child_field: str
    edge_name: str
    edge_direction_to_parent: bool
    remove_parent_value_on_join: bool


@dataclass
class PipelineRecord[T: PyeedBase]:
    data: T
    children: list[ChildRecord[T]] = field(default_factory=list)
    embeddings: dict[str, np.ndarray] = field(default_factory=dict)


@dataclass
class StageConfig:
    """Configuration for a pipeline stage."""

    stage: PipelineStage
    input_queues: list[str]
    output_queues: list[str]
    task_id: TaskID | None


class Pipeline:
    def __init__(
        self,
        progress: Progress | None = None,
    ):
        self.stages: list[StageConfig] = []
        self.queues: dict[str, asyncio.Queue[Any]] = {}
        self.progress = progress or create_progress()
        self.context = PipelineContext(progress=self.progress)

    def add_queue(self, name: str, maxsize: int) -> str:
        """Register a named queue."""
        self.queues[name] = asyncio.Queue(maxsize=maxsize)
        return name

    def add_stage(
        self,
        stage: PipelineStage,
        input_queues: list[str],
        output_queues: list[str],
        task_id: TaskID | None = None,
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
            task_id=task_id,
        )
        self.stages.append(config)

    async def run(
        self,
    ) -> dict[str, int]:
        """Run all stages concurrently."""

        # Launch all stages
        tasks = []
        for config in self.stages:
            in_qs = {name: self.queues[name] for name in config.input_queues}
            out_qs = {name: self.queues[name] for name in config.output_queues}

            tasks.append(
                asyncio.create_task(
                    config.stage.run(
                        input_queues=in_qs,
                        output_queues=out_qs,
                        context=self.context,
                        progress=self.progress,
                        task_id=config.task_id,
                    )
                )
            )

        await asyncio.gather(*tasks)
        return self.context.stats
