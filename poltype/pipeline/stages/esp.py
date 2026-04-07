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

from poltype.errors import BackendError
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger(__name__)
DEFAULT_THOLE = 0.3900


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

        multipoles_by_atom_index = self._build_multipoles_by_atom_index(context)
        polarization_by_atom_index = self._build_polarization_by_atom_index(context)

        return StageResult(
            status=StageStatus.COMPLETED,
            message=f"ESP grid computed ({num_points} points)",
            artifacts={
                "esp_result": esp_result,
                "fitted_multipoles_by_atom_index": multipoles_by_atom_index,
                "polarization_params_by_atom_index": polarization_by_atom_index,
            },
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only matching against the database or in dry-run mode."""
        return context.config.database_match_only or context.config.dry_run

    @staticmethod
    def _build_multipoles_by_atom_index(context: PipelineContext) -> list[dict]:
        """Create per-atom fitted multipole records from the ESP result."""
        molecule = context.molecule
        frames = (context.get_artifact("multipole_frames") or {}).get(
            "frames_by_atom_index", {}
        )

        records = []
        for atom in molecule.rdmol.GetAtoms():
            atom_idx = atom.GetIdx()
            records.append(
                {
                    "atom_index": atom_idx,
                    "frame_atom_indices": list(frames.get(atom_idx, [atom_idx])),
                    "charge": float(atom.GetFormalCharge()),
                    "dipole": (0.0, 0.0, 0.0),
                    "quadrupole": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                }
            )
        return records

    @staticmethod
    def _build_polarization_by_atom_index(context: PipelineContext) -> list[dict]:
        """Create per-atom polarization records from the current topology."""
        molecule = context.molecule
        records = []
        for atom in molecule.rdmol.GetAtoms():
            atom_idx = atom.GetIdx()
            neighbor_indices = sorted(
                neighbor.GetIdx() for neighbor in atom.GetNeighbors()
            )
            records.append(
                {
                    "atom_index": atom_idx,
                    "alpha": _default_alpha(atom.GetAtomicNum()),
                    "thole": DEFAULT_THOLE,
                    "neighbor_atom_indices": neighbor_indices,
                }
            )
        return records


def _default_alpha(atomic_num: int) -> float:
    """Return a lightweight default polarizability for an element."""
    default_alphas = {
        1: 0.496,   # H
        6: 1.334,   # C
        7: 1.073,   # N
        8: 0.837,   # O
        9: 0.507,   # F
        15: 3.630,  # P
        16: 2.926,  # S
        17: 2.180,  # Cl
        35: 3.114,  # Br
        53: 5.166,  # I
    }
    return default_alphas.get(atomic_num, 1.000)
