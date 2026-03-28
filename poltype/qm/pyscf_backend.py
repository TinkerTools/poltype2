"""
poltype.qm.pyscf_backend – PySCF QM backend adapter (stub).

Phase 1 status: **stub adapter** – the class satisfies the
:class:`~poltype.qm.backend.QMBackend` interface.  The existing PySCF
code paths remain active until Phase 4.
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


class PySCFBackend(QMBackend):
    """PySCF-based QM backend.

    Parameters
    ----------
    work_dir:
        Scratch directory for PySCF output files.
    num_proc:
        Number of OpenMP threads (``lib.num_threads``).
    max_mem_gb:
        Memory limit in GB (``mol.max_memory``).
    """

    def __init__(
        self,
        work_dir: Optional[Path] = None,
        num_proc: int = 1,
        max_mem_gb: int = 16,
    ) -> None:
        super().__init__(work_dir=work_dir)
        self.num_proc = num_proc
        self.max_mem_gb = max_mem_gb

    @property
    def name(self) -> str:
        return "pyscf"

    def optimize_geometry(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> OptimizationResult:
        """Geometry optimisation via PySCF.  Phase 1 stub."""
        raise NotImplementedError(
            "PySCFBackend.optimize_geometry is a Phase 1 stub. "
            "Not yet implemented."
        )

    def compute_esp_grid(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "aug-cc-pVTZ",
        **kwargs,
    ) -> ESPGridResult:
        """ESP grid calculation via PySCF.  Phase 1 stub."""
        raise NotImplementedError(
            "PySCFBackend.compute_esp_grid is a Phase 1 stub. "
            "Not yet implemented."
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
        """Relaxed torsion scan via PySCF.  Phase 1 stub."""
        raise NotImplementedError(
            "PySCFBackend.torsion_scan is a Phase 1 stub. "
            "Not yet implemented."
        )

    def compute_wbo_matrix(
        self,
        molecule: Molecule,
        method: str = "HF",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> WBOResult:
        """Wiberg bond order matrix via PySCF.  Phase 1 stub."""
        raise NotImplementedError(
            "PySCFBackend.compute_wbo_matrix is a Phase 1 stub. "
            "Not yet implemented."
        )
