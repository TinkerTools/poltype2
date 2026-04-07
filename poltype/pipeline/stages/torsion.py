"""
poltype.pipeline.stages.torsion – torsion scanning and fitting stage.

Delegates to :meth:`QMBackend.torsion_scan` for each rotatable bond
in the molecule and stores the scan results as pipeline artifacts.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Tuple

import numpy as np

from poltype.errors import BackendError
from poltype.output.params import TorsionFold, TorsionParam
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

        matched_torsions = set()
        db_result = context.get_artifact("db_match_result")
        if db_result is not None:
            matched_torsions = set(db_result.matched_torsions)

        bonds_to_scan = [
            bond for bond in rotatable
            if not self._torsion_covered(bond, matched_torsions)
        ]

        if not bonds_to_scan:
            logger.info("All rotatable bonds were matched from the database")
            return StageResult(
                status=StageStatus.COMPLETED,
                message="All rotatable bonds matched from database",
                artifacts={
                    "torsion_scans": {},
                    "torsion_params": [],
                },
            )

        qm = context.config.qm
        scans: Dict[Tuple[int, int], object] = {}

        for bond in bonds_to_scan:
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

        torsion_params = self._build_torsion_params(context, scans)
        logger.info(
            "Torsion scanning complete: %d bonds scanned", len(scans),
        )

        return StageResult(
            status=StageStatus.COMPLETED,
            message=f"Scanned {len(scans)} rotatable bond(s)",
            artifacts={
                "torsion_scans": scans,
                "torsion_params": torsion_params,
            },
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when torsion fitting is disabled, database-only, or dry-run mode."""
        return (
            not context.config.do_torsion_fit
            or context.config.database_match_only
            or context.config.dry_run
        )

    @staticmethod
    def _torsion_covered(
        bond: tuple[int, int],
        matched_torsions: set[tuple[int, int, int, int]],
    ) -> bool:
        for torsion in matched_torsions:
            if len(torsion) < 4:
                continue
            central = (torsion[1], torsion[2])
            if central == bond or central == (bond[1], bond[0]):
                return True
        return False

    @staticmethod
    def _build_torsion_params(
        context: PipelineContext,
        scans: Dict[Tuple[int, int], object],
    ) -> list[TorsionParam]:
        type_records = context.get_artifact("type_records") or []
        type_map = {
            record.atom_index: (record.atom_type if record.atom_type != 0 else record.atom_index + 1)
            for record in type_records
        }

        params_by_type: dict[tuple[int, int, int, int], TorsionParam] = {}
        rdmol = context.molecule.rdmol
        for bond, scan_result in scans.items():
            torsion_indices = _torsion_indices_for_bond(rdmol, bond)
            if torsion_indices is None:
                continue

            atom_types = tuple(type_map.get(idx, idx + 1) for idx in torsion_indices)
            params_by_type.setdefault(
                atom_types,
                TorsionParam(
                    atom_types=atom_types,
                    folds=_fit_torsion_folds(scan_result),
                ),
            )

        return sorted(params_by_type.values(), key=lambda tp: tp.atom_types)


def _torsion_indices_for_bond(rdmol, bond: tuple[int, int]) -> tuple[int, int, int, int] | None:
    j, k = bond
    left = sorted(
        neighbor.GetIdx()
        for neighbor in rdmol.GetAtomWithIdx(j).GetNeighbors()
        if neighbor.GetIdx() != k
    )
    right = sorted(
        neighbor.GetIdx()
        for neighbor in rdmol.GetAtomWithIdx(k).GetNeighbors()
        if neighbor.GetIdx() != j
    )
    if not left or not right:
        return None
    return (left[0], j, k, right[0])


def _fit_torsion_folds(scan_result) -> list[TorsionFold]:
    energies = np.asarray(scan_result.energies, dtype=float)
    if energies.size == 0:
        amplitude = 0.0
    else:
        amplitude = (energies.max() - energies.min()) * 627.509474 / 2.0
    return [TorsionFold(amplitude=float(amplitude), phase=0.0, periodicity=3)]
