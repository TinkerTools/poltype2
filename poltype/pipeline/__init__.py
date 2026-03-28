"""
poltype.pipeline – composable pipeline orchestration layer.

Re-exports the core framework types so callers can write::

    from poltype.pipeline import Stage, PipelineContext, PipelineRunner
"""

from __future__ import annotations

from poltype.pipeline.checkpoint import CheckpointData, CheckpointManager
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.events import EventBus, EventData, PipelineEvent
from poltype.pipeline.factory import build_default_pipeline
from poltype.pipeline.runner import PipelineRunner
from poltype.pipeline.stage import Stage, StageResult, StageStatus

__all__ = [
    "Stage",
    "StageResult",
    "StageStatus",
    "PipelineContext",
    "PipelineRunner",
    "EventBus",
    "EventData",
    "PipelineEvent",
    "CheckpointData",
    "CheckpointManager",
    "build_default_pipeline",
]
