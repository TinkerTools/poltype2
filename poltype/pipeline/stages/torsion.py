"""
poltype.pipeline.stages.torsion – torsion scanning and fitting stage.

Stub that will delegate to :meth:`QMBackend.torsion_scan` and torsion
fitting routines in Phase 4.
"""

from __future__ import annotations

from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus


class TorsionFittingStage(Stage):
    """Scan and fit torsion parameters.

    .. note::
        This is a **stub** — the real implementation will be added in
        Phase 4.
    """

    def __init__(self) -> None:
        super().__init__(name="torsion_fitting")

    def execute(self, context: PipelineContext) -> StageResult:
        """Stub: return COMPLETED with an informational message."""
        return StageResult(
            status=StageStatus.COMPLETED,
            message=(
                "stub — will delegate to QMBackend.torsion_scan() "
                "and torsion fitting in Phase 4"
            ),
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when torsion fitting is disabled or database-only mode."""
        return (
            not context.config.do_torsion_fit
            or context.config.database_match_only
        )
