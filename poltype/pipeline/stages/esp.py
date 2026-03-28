"""
poltype.pipeline.stages.esp – ESP grid computation and fitting stage.

Delegates to :meth:`QMBackend.compute_esp_grid` to compute the
electrostatic potential on a grid around the molecule and stores the
result as a pipeline artifact.
"""

from __future__ import annotations

import logging

from poltype.errors import BackendError, StageError
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger(__name__)


class ESPFittingStage(Stage):
    """Compute the ESP grid and fit electrostatic parameters.

    Calls :meth:`QMBackend.compute_esp_grid` with the ESP method and
    basis set from ``context.config.qm``.
    """

    def __init__(self) -> None:
        super().__init__(name="esp_fitting")

    def execute(self, context: PipelineContext) -> StageResult:
        """Compute the ESP grid via the configured QM backend.

        Returns
        -------
        StageResult
            ``COMPLETED`` with an ``esp_result`` artifact on success,
            ``FAILED`` on error.
        """
        if context.molecule is None:
            return StageResult(
                status=StageStatus.FAILED,
                message="No molecule in context",
            )

        if context.backend is None:
            return StageResult(
                status=StageStatus.FAILED,
                message="No QM backend in context",
            )

        qm = context.config.qm
        molecule = context.molecule

        logger.info(
            "Computing ESP grid: method=%s basis=%s",
            qm.esp_method,
            qm.esp_basis_set,
        )

        try:
            esp_result = context.backend.compute_esp_grid(
                molecule,
                method=qm.esp_method,
                basis_set=qm.esp_basis_set,
            )
        except NotImplementedError:
            raise
        except Exception as exc:
            raise BackendError(
                f"ESP grid computation failed: {exc}",
                backend_name=context.backend.name,
            ) from exc

        num_points = len(esp_result.esp_values)
        logger.info("ESP grid computed: %d grid points", num_points)

        return StageResult(
            status=StageStatus.COMPLETED,
            message=f"ESP grid computed ({num_points} points)",
            artifacts={"esp_result": esp_result},
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only matching against the database or in dry-run mode."""
        return context.config.database_match_only or context.config.dry_run
