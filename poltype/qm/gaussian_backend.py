"""
poltype.qm.gaussian_backend – Gaussian QM backend adapter (stub).

Phase 1 status: **stub adapter** – the class satisfies the
:class:`~poltype.qm.backend.QMBackend` interface.  The legacy
``use_gaus`` / ``use_gausoptonly`` code paths in
``PoltypeModules/optimization.py`` and related modules remain the
active implementation until Phase 4.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import numpy as np

from poltype.molecule.molecule import Molecule
from poltype.qm.backend import (
    ESPGridResult,
    OptimizationResult,
    QMBackend,
    TorsionScanResult,
    WBOResult,
)


class GaussianBackend(QMBackend):
    """Gaussian-based QM backend.

    Parameters
    ----------
    work_dir:
        Scratch directory for Gaussian input/output files.
    gaus_exe:
        Path / name of the Gaussian executable (e.g. ``"g16"``).
    num_proc:
        Number of parallel threads (``%NProcShared`` in the input).
    max_mem_gb:
        Memory limit in GB (``%Mem`` in the input).
    """

    def __init__(
        self,
        work_dir: Optional[Path] = None,
        gaus_exe: str = "g16",
        num_proc: int = 1,
        max_mem_gb: int = 16,
    ) -> None:
        super().__init__(work_dir=work_dir)
        self.gaus_exe = gaus_exe
        self.num_proc = num_proc
        self.max_mem_gb = max_mem_gb

    @property
    def name(self) -> str:
        return "gaussian"

    def optimize_geometry(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> OptimizationResult:
        """Geometry optimisation via Gaussian.  Phase 1 stub."""
        raise NotImplementedError(
            "GaussianBackend.optimize_geometry is a Phase 1 stub. "
            "Use PoltypeModules/optimization.py (use_gaus=True) directly."
        )

    def compute_esp_grid(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "aug-cc-pVTZ",
        **kwargs,
    ) -> ESPGridResult:
        """ESP grid calculation via Gaussian.  Phase 1 stub."""
        raise NotImplementedError(
            "GaussianBackend.compute_esp_grid is a Phase 1 stub. "
            "Use PoltypeModules/electrostaticpotential.py directly."
        )

    def torsion_scan(
        self,
        molecule: Molecule,
        rotatable_bond: Tuple[int, int],
        method: str = "wB97X-D",
        basis_set: str = "6-31G*",
        angles: Optional[np.ndarray] = None,
        **kwargs,
    ) -> TorsionScanResult:
        """Relaxed torsion scan via Gaussian.  Phase 1 stub."""
        raise NotImplementedError(
            "GaussianBackend.torsion_scan is a Phase 1 stub. "
            "Use PoltypeModules/torsiongenerator.py directly."
        )

    def compute_wbo_matrix(
        self,
        molecule: Molecule,
        method: str = "HF",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> WBOResult:
        """Wiberg bond order matrix via Gaussian.  Phase 1 stub."""
        raise NotImplementedError(
            "GaussianBackend.compute_wbo_matrix is a Phase 1 stub. "
            "Use PoltypeModules/fragmenter.py directly."
        )
