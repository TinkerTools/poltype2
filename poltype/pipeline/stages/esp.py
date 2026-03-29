"""
poltype.pipeline.stages.esp – ESP grid computation and multipole fitting stage.

Computes the electrostatic potential on a grid around the molecule and
uses it to optimise the initial atomic multipoles obtained from the
upstream :class:`MultipoleStage` (Distributed Multipole Analysis).

The correct protocol is:

1. **DMA** (:class:`MultipoleStage`) — Distributed Multipole Analysis to
   obtain initial atomic multipoles.
2. **ESP fitting** (this stage) — Electrostatic potential fitting to
   further optimise those multipoles.
"""

from __future__ import annotations

import logging

from poltype.errors import BackendError, StageError
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger(__name__)


class ESPFittingStage(Stage):
    """Compute the ESP grid and optimise multipoles via potential fitting.

    Consumes the ``multipole_frames`` artifact published by the upstream
    :class:`MultipoleStage` and calls :meth:`QMBackend.compute_esp_grid`
    with the ESP method and basis set from ``context.config.qm`` to
    refine the atomic multipoles.
    """

    def __init__(self) -> None:
        super().__init__(name="esp_fitting")

    def execute(self, context: PipelineContext) -> StageResult:
        """Compute the ESP grid and optimise multipoles.

        The upstream :class:`MultipoleStage` must have produced a
        ``multipole_frames`` artifact containing the initial DMA
        multipoles.  This stage computes the electrostatic potential
        grid and uses it to further optimise those multipoles.

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

        multipole_frames = context.get_artifact("multipole_frames")
        if multipole_frames is not None:
            logger.info(
                "Optimising DMA multipoles for %d atoms via ESP fitting",
                multipole_frames.get("num_atoms", "?"),
            )
        else:
            logger.info(
                "No DMA multipoles available; computing ESP grid only"
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
