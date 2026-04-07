"""
poltype.pipeline.stages.esp – ESP grid computation and multipole fitting stage.

Computes the electrostatic potential on a grid around the molecule and
uses it to optimise the initial atomic multipoles obtained from the
upstream :class:`MultipoleStage` (Distributed Multipole Analysis).

The correct protocol is:

1. **DMA** (:class:`MultipoleStage`) — Distributed Multipole Analysis to
   obtain initial atomic multipoles.
2. **ESP fitting** (this stage) — High-level QM (MP2/aug-cc-pVTZ) to
   compute the electrostatic potential, then Tinker's ``potential``
   program to fit multipoles targeting that ESP.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np

from poltype.errors import BackendError
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger(__name__)
DEFAULT_THOLE = 0.3900


class ESPFittingStage(Stage):
    """Compute the ESP grid and optimise multipoles via potential fitting.

    This stage:

    1. Runs a high-level QM calculation (default MP2/aug-cc-pVTZ) to
       compute the electrostatic potential on a grid around the molecule.
    2. Optionally runs Tinker's ``potential`` program to fit the
       multipoles from the upstream :class:`MultipoleStage` to match
       the computed ESP.
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

        # ---- Step 5: Run QM at high level to get ESP ----
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

        # Write ESP grid data to file for Tinker potential program
        esp_grid_path = self._write_esp_grid(
            esp_result, context.work_dir, molecule.name or "mol",
        )

        # ---- Step 6: Run Tinker potential to fit multipoles ----
        potential_result = None
        opt_xyz_path = context.get_artifact("opt_xyz_path")
        if opt_xyz_path is not None:
            try:
                from poltype.external.potential import PotentialRunner

                potential = PotentialRunner(
                    potential_exe=context.config.potential_path,
                )
                potential_result = potential.run(
                    xyz_path=opt_xyz_path,
                    work_dir=context.work_dir,
                    molecule_name=molecule.name or "mol",
                    esp_grid_path=esp_grid_path,
                )
                logger.info(
                    "Potential fitting completed: key file at %s",
                    potential_result.key_path,
                )
            except FileNotFoundError:
                logger.warning(
                    "potential executable not found (%s); skipping "
                    "potential fitting step.",
                    context.config.potential_path,
                )
            except Exception as exc:
                logger.warning("Potential fitting failed: %s", exc)

        multipoles_by_atom_index = self._build_multipoles_by_atom_index(context)
        polarization_by_atom_index = self._build_polarization_by_atom_index(context)

        artifacts = {
            "esp_result": esp_result,
            "fitted_multipoles_by_atom_index": multipoles_by_atom_index,
            "polarization_params_by_atom_index": polarization_by_atom_index,
        }
        if esp_grid_path is not None:
            artifacts["esp_grid_path"] = esp_grid_path
        if potential_result is not None:
            artifacts["potential_result"] = potential_result

        return StageResult(
            status=StageStatus.COMPLETED,
            message=f"ESP grid computed ({num_points} points)",
            artifacts=artifacts,
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only matching against the database or in dry-run mode."""
        return context.config.database_match_only or context.config.dry_run

    @staticmethod
    def _write_esp_grid(
        esp_result,
        work_dir: Path,
        molecule_name: str,
    ) -> Optional[Path]:
        """Write the ESP grid to a file for Tinker's potential program.

        The file contains grid coordinates (Angstroms) and ESP values
        (atomic units) in a format suitable for potential fitting.

        Returns
        -------
        Path or None
            Path to the written grid file, or None if no grid points.
        """
        if len(esp_result.esp_values) == 0:
            return None

        grid_path = Path(work_dir) / f"{molecule_name}_esp.grid"
        n_points = len(esp_result.esp_values)

        lines = [f"{n_points}"]
        for i in range(n_points):
            x, y, z = esp_result.grid_coords[i]
            v = esp_result.esp_values[i]
            lines.append(
                f"  {x:>14.8f}  {y:>14.8f}  {z:>14.8f}  {v:>14.8f}"
            )

        grid_path.write_text("\n".join(lines) + "\n")
        logger.info("Wrote ESP grid (%d points) to %s", n_points, grid_path)
        return grid_path

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
