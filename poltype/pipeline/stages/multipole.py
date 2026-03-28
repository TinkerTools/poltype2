"""
poltype.pipeline.stages.multipole – multipole / GDMA stage.

Stub that will delegate to GDMA / poledit in Phase 4.
"""

from __future__ import annotations

from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus


class MultipoleStage(Stage):
    """Compute distributed multipoles via GDMA and poledit.

    .. note::
        This is a **stub** — the real implementation will be added in
        Phase 4.
    """

    def __init__(self) -> None:
        super().__init__(name="multipole")

    def execute(self, context: PipelineContext) -> StageResult:
        """Stub: return COMPLETED with an informational message."""
        return StageResult(
            status=StageStatus.COMPLETED,
            message="stub — will delegate to GDMA/poledit in Phase 4",
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only matching against the database."""
        return context.config.database_match_only
