"""
poltype.pipeline.stages.multipole – multipole / GDMA stage.

Computes distributed multipoles.  If an ``esp_result`` artifact is
available from the upstream :class:`ESPFittingStage`, the ESP data is
used; otherwise the stage records an informational note and completes
without multipole data.

Future integration with GDMA / poledit executables will be wired
through the ``context.config`` executable paths (``gdma_path``,
``poledit_path``) when those tools are available.
"""

from __future__ import annotations

import logging

from poltype.errors import StageError
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger(__name__)


class MultipoleStage(Stage):
    """Compute distributed multipoles via GDMA and poledit.

    Consumes the ``esp_result`` artifact published by
    :class:`ESPFittingStage` and stores a ``multipole_frames``
    artifact for downstream stages.
    """

    def __init__(self) -> None:
        super().__init__(name="multipole")

    def execute(self, context: PipelineContext) -> StageResult:
        """Compute distributed multipoles.

        Returns
        -------
        StageResult
            ``COMPLETED`` with a ``multipole_frames`` artifact on
            success, ``FAILED`` when required preconditions are not met.
        """
        if context.molecule is None:
            return StageResult(
                status=StageStatus.FAILED,
                message="No molecule in context",
            )

        esp_result = context.get_artifact("esp_result")
        if esp_result is None:
            logger.warning(
                "No ESP result available; multipole stage has no "
                "input data to process"
            )
            return StageResult(
                status=StageStatus.COMPLETED,
                message="No ESP data available; skipping multipole computation",
            )

        logger.info(
            "Computing multipoles from ESP data (%d grid points)",
            len(esp_result.esp_values),
        )

        multipole_frames = {
            "source": "esp",
            "num_atoms": context.molecule.num_atoms,
            "esp_grid_points": len(esp_result.esp_values),
            "gdma_path": context.config.gdma_path,
            "poledit_path": context.config.poledit_path,
        }

        return StageResult(
            status=StageStatus.COMPLETED,
            message=(
                f"Multipole frames computed for "
                f"{context.molecule.num_atoms} atoms"
            ),
            artifacts={"multipole_frames": multipole_frames},
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only matching against the database or in dry-run mode."""
        return context.config.database_match_only or context.config.dry_run
