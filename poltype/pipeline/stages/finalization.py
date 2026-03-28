"""
poltype.pipeline.stages.finalization – final output assembly stage.

Stub that will copy the final ``.xyz`` / ``.key`` files in Phase 4.
"""

from __future__ import annotations

from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus


class FinalizationStage(Stage):
    """Assemble final output files (.xyz, .key).

    .. note::
        This is a **stub** — the real implementation will be added in
        Phase 4.
    """

    def __init__(self) -> None:
        super().__init__(name="finalization")

    def execute(self, context: PipelineContext) -> StageResult:
        """Stub: return COMPLETED with an informational message."""
        return StageResult(
            status=StageStatus.COMPLETED,
            message="stub — will copy final .xyz/.key files in Phase 4",
        )
