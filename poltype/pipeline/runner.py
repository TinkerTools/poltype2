"""
poltype.pipeline.runner – pipeline orchestrator.

The :class:`PipelineRunner` iterates an ordered list of :class:`Stage`
objects, executing each in turn against a shared
:class:`PipelineContext`.  Stages may be skipped (config flags),
succeed, or fail — a failure halts the pipeline immediately.

This module has **no dependency on PoltypeModules** and can be
unit-tested in isolation.

Example::

    from poltype.pipeline import PipelineRunner, PipelineContext
    from poltype.pipeline.stages import InputPreparationStage

    runner = PipelineRunner()
    runner.add_stage(InputPreparationStage())
    ctx = runner.run(PipelineContext(config=my_config))
"""

from __future__ import annotations

import time
from typing import List, Optional

from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus


class PipelineRunner:
    """Executes a sequence of pipeline stages against a shared context.

    Parameters
    ----------
    stages:
        Optional initial list of stages.  More stages can be added
        later via :meth:`add_stage`.
    """

    def __init__(self, stages: Optional[List[Stage]] = None) -> None:
        self._stages: List[Stage] = list(stages) if stages else []

    def add_stage(self, stage: Stage) -> PipelineRunner:
        """Append a stage to the pipeline.

        Returns ``self`` so calls can be chained::

            runner.add_stage(A()).add_stage(B())
        """
        self._stages.append(stage)
        return self

    def run(self, context: PipelineContext) -> PipelineContext:
        """Execute all stages in order.

        For each stage the runner:

        1. Checks :meth:`Stage.should_skip` — if ``True``, records a
           ``SKIPPED`` result and moves on.
        2. Calls :meth:`Stage.execute` and records the result.
        3. If the result status is ``FAILED``, stops the pipeline.
        4. Merges ``result.artifacts`` into ``context.artifacts``.

        Parameters
        ----------
        context:
            Mutable pipeline state shared across stages.

        Returns
        -------
        PipelineContext
            The same context object, enriched with stage results and
            artifacts.
        """
        for stage in self._stages:
            # -- skip check --
            if stage.should_skip(context):
                result = StageResult(
                    status=StageStatus.SKIPPED,
                    message=f"Stage '{stage.name}' skipped by should_skip()",
                )
                context.stage_results[stage.name] = result
                continue

            # -- execute --
            t0 = time.monotonic()
            result = stage.execute(context)
            result.elapsed_seconds = time.monotonic() - t0

            context.stage_results[stage.name] = result

            # -- stop on failure --
            if result.status is StageStatus.FAILED:
                break

            # -- merge artifacts --
            context.artifacts.update(result.artifacts)

        return context

    def __repr__(self) -> str:
        return f"PipelineRunner(stages={len(self._stages)})"
