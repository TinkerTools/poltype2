"""
poltype.pipeline.stages.geometry_opt – QM geometry optimisation stage.

Delegates to :meth:`QMBackend.optimize_geometry` to run a QM geometry
optimisation and updates the molecule in the pipeline context with the
optimised coordinates.
"""

from __future__ import annotations

import logging

from poltype.errors import BackendError, StageError
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger(__name__)


class GeometryOptimizationStage(Stage):
    """Run QM geometry optimisation on the current molecule.

    Calls :meth:`QMBackend.optimize_geometry` with the method and
    basis set from ``context.config.qm``, then replaces the molecule
    in the context with one carrying the optimised coordinates.
    """

    def __init__(self) -> None:
        super().__init__(name="geometry_optimization")

    def execute(self, context: PipelineContext) -> StageResult:
        """Optimise geometry via the configured QM backend.

        Returns
        -------
        StageResult
            ``COMPLETED`` with ``opt_result`` and ``molecule`` artifacts
            on success, ``FAILED`` on error.
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
            "Running geometry optimisation: method=%s basis=%s",
            qm.opt_method,
            qm.opt_basis_set,
        )

        try:
            opt_result = context.backend.optimize_geometry(
                molecule,
                method=qm.opt_method,
                basis_set=qm.opt_basis_set,
            )
        except NotImplementedError:
            raise
        except Exception as exc:
            raise BackendError(
                f"Geometry optimisation failed: {exc}",
                backend_name=context.backend.name,
            ) from exc

        if not opt_result.converged:
            return StageResult(
                status=StageStatus.FAILED,
                message="Geometry optimisation did not converge",
                artifacts={"opt_result": opt_result},
            )

        new_mol = molecule.with_updated_coordinates(opt_result.coordinates)
        context.update_molecule(new_mol)

        logger.info(
            "Geometry optimisation converged: energy=%.6f Hartree",
            opt_result.energy,
        )

        return StageResult(
            status=StageStatus.COMPLETED,
            message=f"Optimised geometry (energy={opt_result.energy:.6f} Ha)",
            artifacts={
                "opt_result": opt_result,
                "molecule": new_mol,
            },
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only matching against the database or in dry-run mode."""
        return context.config.database_match_only or context.config.dry_run
