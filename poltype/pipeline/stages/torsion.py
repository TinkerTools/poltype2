"""
poltype.pipeline.stages.torsion – torsion scanning and fitting stage.

Delegates to :meth:`QMBackend.torsion_scan` for each rotatable bond
in the molecule and stores the scan results as pipeline artifacts.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Tuple

from poltype.errors import BackendError, StageError
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger(__name__)


class TorsionFittingStage(Stage):
    """Scan and fit torsion parameters.

    Iterates over the rotatable bonds of the current molecule, runs
    :meth:`QMBackend.torsion_scan` for each, and stores the collected
    results as the ``torsion_scans`` artifact.
    """

    def __init__(self) -> None:
        super().__init__(name="torsion_fitting")

    def execute(self, context: PipelineContext) -> StageResult:
        """Run torsion scans for all rotatable bonds.

        Returns
        -------
        StageResult
            ``COMPLETED`` with a ``torsion_scans`` artifact on success,
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

        molecule = context.molecule
        rotatable = molecule.rotatable_bonds

        if not rotatable:
            logger.info("No rotatable bonds found; nothing to scan")
            return StageResult(
                status=StageStatus.COMPLETED,
                message="No rotatable bonds to scan",
                artifacts={"torsion_scans": {}},
            )

        qm = context.config.qm
        scans: Dict[Tuple[int, int], object] = {}

        for bond in rotatable:
            logger.info(
                "Scanning torsion around bond %s: method=%s basis=%s",
                bond,
                qm.tor_opt_method,
                qm.tor_opt_basis_set,
            )
            try:
                scan_result = context.backend.torsion_scan(
                    molecule,
                    rotatable_bond=bond,
                    method=qm.tor_opt_method,
                    basis_set=qm.tor_opt_basis_set,
                )
            except NotImplementedError:
                raise
            except Exception as exc:
                raise BackendError(
                    f"Torsion scan failed for bond {bond}: {exc}",
                    backend_name=context.backend.name,
                ) from exc

            scans[bond] = scan_result

        logger.info(
            "Torsion scanning complete: %d bonds scanned", len(scans),
        )

        return StageResult(
            status=StageStatus.COMPLETED,
            message=f"Scanned {len(scans)} rotatable bond(s)",
            artifacts={"torsion_scans": scans},
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when torsion fitting is disabled, database-only, or dry-run mode."""
        return (
            not context.config.do_torsion_fit
            or context.config.database_match_only
            or context.config.dry_run
        )
