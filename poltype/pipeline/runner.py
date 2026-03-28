"""
poltype.pipeline.runner – pipeline orchestrator.

The :class:`PipelineRunner` iterates an ordered list of :class:`Stage`
objects, executing each in turn against a shared
:class:`PipelineContext`.  Stages may be skipped (config flags),
succeed, or fail — a failure halts the pipeline immediately.

Optionally, an :class:`EventBus` can be attached to broadcast
:class:`PipelineEvent` notifications, and a :class:`CheckpointManager`
can be provided to persist progress and enable resume.

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

import logging
import time
from typing import List, Optional

from poltype.pipeline.checkpoint import CheckpointManager
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.events import EventBus, EventData, PipelineEvent
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger("poltype.pipeline.runner")


class PipelineRunner:
    """Executes a sequence of pipeline stages against a shared context.

    Parameters
    ----------
    stages:
        Optional initial list of stages.  More stages can be added
        later via :meth:`add_stage`.
    event_bus:
        Optional :class:`EventBus` for broadcasting pipeline events.
    checkpoint_manager:
        Optional :class:`CheckpointManager` for persisting progress
        and enabling resume.
    """

    def __init__(
        self,
        stages: Optional[List[Stage]] = None,
        event_bus: Optional[EventBus] = None,
        checkpoint_manager: Optional[CheckpointManager] = None,
    ) -> None:
        self._stages: List[Stage] = list(stages) if stages else []
        self._event_bus: Optional[EventBus] = event_bus
        self._checkpoint: Optional[CheckpointManager] = checkpoint_manager

    @property
    def event_bus(self) -> Optional[EventBus]:
        """The attached event bus (if any)."""
        return self._event_bus

    @event_bus.setter
    def event_bus(self, value: Optional[EventBus]) -> None:
        self._event_bus = value

    @property
    def checkpoint_manager(self) -> Optional[CheckpointManager]:
        """The attached checkpoint manager (if any)."""
        return self._checkpoint

    @checkpoint_manager.setter
    def checkpoint_manager(self, value: Optional[CheckpointManager]) -> None:
        self._checkpoint = value

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

        1. Checks the checkpoint — if the stage already completed in a
           previous run, records a ``SKIPPED`` result and moves on.
        2. Checks :meth:`Stage.should_skip` — if ``True``, records a
           ``SKIPPED`` result and moves on.
        3. Emits a ``STAGE_STARTED`` event.
        4. Calls :meth:`Stage.execute` and records the result.
        5. If the result status is ``FAILED``, emits ``STAGE_FAILED``,
           saves the checkpoint, and stops the pipeline.
        6. Merges ``result.artifacts`` into ``context.artifacts``,
           emits ``STAGE_COMPLETED``, and saves the checkpoint.

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
        self._emit(EventData(
            event=PipelineEvent.PIPELINE_STARTED,
            message=f"Pipeline starting with {len(self._stages)} stages",
        ))

        for stage in self._stages:
            # -- checkpoint resume check --
            if self._checkpoint and self._checkpoint.should_skip_stage(stage.name):
                result = StageResult(
                    status=StageStatus.SKIPPED,
                    message=(
                        f"Stage '{stage.name}' skipped (completed in "
                        f"previous run)"
                    ),
                )
                context.stage_results[stage.name] = result
                logger.info(
                    "Stage '%s' skipped (checkpoint resume)", stage.name,
                )
                self._emit(EventData(
                    event=PipelineEvent.STAGE_SKIPPED,
                    stage_name=stage.name,
                    message=result.message,
                ))
                continue

            # -- config skip check --
            if stage.should_skip(context):
                result = StageResult(
                    status=StageStatus.SKIPPED,
                    message=f"Stage '{stage.name}' skipped by should_skip()",
                )
                context.stage_results[stage.name] = result
                logger.info(
                    "Stage '%s' skipped", stage.name,
                )
                self._emit(EventData(
                    event=PipelineEvent.STAGE_SKIPPED,
                    stage_name=stage.name,
                    message=result.message,
                ))
                if self._checkpoint:
                    self._checkpoint.save(stage.name, result)
                continue

            # -- execute --
            self._emit(EventData(
                event=PipelineEvent.STAGE_STARTED,
                stage_name=stage.name,
                message=f"Stage '{stage.name}' starting",
            ))
            logger.info("Stage '%s' starting", stage.name)
            t0 = time.monotonic()
            result = stage.execute(context)
            result.elapsed_seconds = time.monotonic() - t0

            context.stage_results[stage.name] = result

            # -- stop on failure --
            if result.status is StageStatus.FAILED:
                logger.error(
                    "Stage '%s' failed (%.2fs): %s",
                    stage.name,
                    result.elapsed_seconds,
                    result.message,
                )
                self._emit(EventData(
                    event=PipelineEvent.STAGE_FAILED,
                    stage_name=stage.name,
                    message=result.message,
                    extra={"elapsed_seconds": result.elapsed_seconds},
                ))
                if self._checkpoint:
                    self._checkpoint.save(stage.name, result)
                break

            # -- merge artifacts --
            context.artifacts.update(result.artifacts)
            logger.info(
                "Stage '%s' completed (%.2fs)",
                stage.name,
                result.elapsed_seconds,
            )
            self._emit(EventData(
                event=PipelineEvent.STAGE_COMPLETED,
                stage_name=stage.name,
                message=result.message,
                extra={"elapsed_seconds": result.elapsed_seconds},
            ))
            if self._checkpoint:
                self._checkpoint.save(stage.name, result)

        self._emit(EventData(
            event=PipelineEvent.PIPELINE_COMPLETED,
            message="Pipeline finished",
        ))

        return context

    # -- internal helpers ----------------------------------------------------

    def _emit(self, data: EventData) -> None:
        """Broadcast an event if an event bus is attached."""
        if self._event_bus is not None:
            self._event_bus.emit(data)

    def __repr__(self) -> str:
        return f"PipelineRunner(stages={len(self._stages)})"
