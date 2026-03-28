"""
poltype.pipeline.stages.geometry_opt – QM geometry optimisation stage.

Stub that will delegate to :meth:`QMBackend.optimize_geometry` in
Phase 4.
"""

from __future__ import annotations

from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus


class GeometryOptimizationStage(Stage):
    """Run QM geometry optimisation on the current molecule.

    .. note::
        This is a **stub** — the real implementation will be added in
        Phase 4.
    """

    def __init__(self) -> None:
        super().__init__(name="geometry_optimization")

    def execute(self, context: PipelineContext) -> StageResult:
        """Stub: return COMPLETED with an informational message."""
        return StageResult(
            status=StageStatus.COMPLETED,
            message=(
                "stub — will delegate to "
                "QMBackend.optimize_geometry() in Phase 4"
            ),
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only matching against the database."""
        return context.config.database_match_only
