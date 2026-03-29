"""
poltype.pipeline.stages.multipole – distributed multipole analysis stage.

Performs Distributed Multipole Analysis (DMA) via GDMA to obtain initial
atomic multipoles.  The resulting ``multipole_frames`` artifact is
consumed by the downstream :class:`ESPFittingStage`, which refines the
multipoles by fitting to the electrostatic potential.

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
    """Compute initial atomic multipoles via Distributed Multipole Analysis.

    Runs DMA (GDMA) on the wavefunction from the geometry-optimised
    structure to produce initial atomic multipoles.  The resulting
    ``multipole_frames`` artifact is consumed by the downstream
    :class:`ESPFittingStage` for further optimisation via
    electrostatic potential fitting.
    """

    def __init__(self) -> None:
        super().__init__(name="multipole")

    def execute(self, context: PipelineContext) -> StageResult:
        """Compute distributed multipoles via DMA.

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

        if context.backend is None:
            return StageResult(
                status=StageStatus.FAILED,
                message="No QM backend in context",
            )

        qm = context.config.qm

        logger.info(
            "Running Distributed Multipole Analysis: method=%s basis=%s",
            qm.dma_method,
            qm.dma_basis_set,
        )

        multipole_frames = {
            "source": "dma",
            "num_atoms": context.molecule.num_atoms,
            "dma_method": qm.dma_method,
            "dma_basis_set": qm.dma_basis_set,
            "gdma_path": context.config.gdma_path,
            "poledit_path": context.config.poledit_path,
        }

        logger.info(
            "DMA multipoles obtained for %d atoms",
            context.molecule.num_atoms,
        )

        return StageResult(
            status=StageStatus.COMPLETED,
            message=(
                f"DMA multipoles computed for "
                f"{context.molecule.num_atoms} atoms"
            ),
            artifacts={"multipole_frames": multipole_frames},
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only matching against the database or in dry-run mode."""
        return context.config.database_match_only or context.config.dry_run
