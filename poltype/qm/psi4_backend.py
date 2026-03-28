"""
poltype.qm.psi4_backend – Psi4 QM backend adapter (stub).

This module provides the :class:`Psi4Backend` class that wraps the
existing Psi4 input-file generation and execution logic found in
``PoltypeModules/optimization.py``, ``electrostaticpotential.py``, and
``torsiongenerator.py``.

Phase 1 status: **stub adapter** – the class satisfies the
:class:`~poltype.qm.backend.QMBackend` interface and delegates to the
legacy module functions.  In later phases the delegation will be
replaced with direct, self-contained Psi4 logic.
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


class Psi4Backend(QMBackend):
    """Psi4-based QM backend.

    Delegates to the legacy ``PoltypeModules/optimization.py`` functions
    in Phase 1.  Subsequent phases will move the Psi4-specific logic
    into this class directly.

    Parameters
    ----------
    work_dir:
        Scratch directory for Psi4 input/output files.
    psi4_args:
        Extra command-line arguments passed to the ``psi4`` executable
        (e.g. ``"--loglevel 30"``).
    num_proc:
        Number of parallel threads for Psi4.
    max_mem_gb:
        Memory limit in GB passed to Psi4 (``memory`` keyword).
    """

    def __init__(
        self,
        work_dir: Optional[Path] = None,
        psi4_args: str = " --loglevel 30 ",
        num_proc: int = 1,
        max_mem_gb: int = 16,
    ) -> None:
        super().__init__(work_dir=work_dir)
        self.psi4_args = psi4_args
        self.num_proc = num_proc
        self.max_mem_gb = max_mem_gb

    @property
    def name(self) -> str:
        return "psi4"

    def optimize_geometry(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> OptimizationResult:
        """Geometry optimisation via Psi4.

        Phase 1 stub: raises ``NotImplementedError`` to signal that the
        delegation to the legacy module is not yet wired up.  Pipeline
        code should continue to call the legacy functions directly until
        Phase 4 replaces this stub with a real implementation.
        """
        raise NotImplementedError(
            "Psi4Backend.optimize_geometry is a Phase 1 stub. "
            "Use PoltypeModules/optimization.py directly for now."
        )

    def compute_esp_grid(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "aug-cc-pVTZ",
        **kwargs,
    ) -> ESPGridResult:
        """ESP grid calculation via Psi4.

        Phase 1 stub.
        """
        raise NotImplementedError(
            "Psi4Backend.compute_esp_grid is a Phase 1 stub. "
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
        """Relaxed torsion scan via Psi4.

        Phase 1 stub.
        """
        raise NotImplementedError(
            "Psi4Backend.torsion_scan is a Phase 1 stub. "
            "Use PoltypeModules/torsiongenerator.py directly."
        )

    def compute_wbo_matrix(
        self,
        molecule: Molecule,
        method: str = "HF",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> WBOResult:
        """Wiberg bond order matrix via Psi4.

        Phase 1 stub.
        """
        raise NotImplementedError(
            "Psi4Backend.compute_wbo_matrix is a Phase 1 stub. "
            "Use PoltypeModules/fragmenter.py directly."
        )
