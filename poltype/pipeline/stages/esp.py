"""
poltype.pipeline.stages.esp – ESP grid computation and fitting stage.

Stub that will delegate to :meth:`QMBackend.compute_esp_grid` in
Phase 4.
"""

from __future__ import annotations

from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus


class ESPFittingStage(Stage):
    """Compute the ESP grid and fit electrostatic parameters.

    .. note::
        This is a **stub** — the real implementation will be added in
        Phase 4.
    """

    def __init__(self) -> None:
        super().__init__(name="esp_fitting")

    def execute(self, context: PipelineContext) -> StageResult:
        """Stub: return COMPLETED with an informational message."""
        return StageResult(
            status=StageStatus.COMPLETED,
            message=(
                "stub — will delegate to "
                "QMBackend.compute_esp_grid() in Phase 4"
            ),
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only matching against the database."""
        return context.config.database_match_only
