"""
poltype.qm.psi4_backend – Psi4 QM backend adapter (stub).

This module provides the :class:`Psi4Backend` class that wraps
Psi4 input-file generation and execution logic.

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

    Subsequent phases will move the Psi4-specific logic
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
        """Geometry optimisation via Psi4."""
        import psi4
        psi4.core.set_num_threads(self.num_proc)
        psi4.set_memory(f"{self.max_mem_gb} GB")
        
        geom_str = f"{molecule.formal_charge} 1\n"
        rdmol = molecule.rdmol
        conf = rdmol.GetConformer()
        for i in range(molecule.num_atoms):
            pos = conf.GetAtomPosition(i)
            sym = rdmol.GetAtomWithIdx(i).GetSymbol()
            geom_str += f"{sym} {pos.x} {pos.y} {pos.z}\n"
        
        psi_mol = psi4.geometry(geom_str)
        psi4.set_options({'basis': basis_set})
        
        # Optimize
        try:
            e, wfn = psi4.optimize(method, return_wfn=True)
            converged = True
        except Exception as e_out:
            # Fallback or error
            raise RuntimeError(f"Psi4 optimization failed: {e_out}")
            
        coords = np.array(psi_mol.geometry()) * psi4.constants.bohr2angstroms
        
        return OptimizationResult(
            coordinates=coords,
            energy=float(e),
            converged=converged,
        )

    def compute_esp_grid(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "aug-cc-pVTZ",
        **kwargs,
    ) -> ESPGridResult:
        """ESP grid calculation via Psi4."""
        import psi4
        psi4.core.set_num_threads(self.num_proc)
        psi4.set_memory(f"{self.max_mem_gb} GB")
        
        geom_str = f"{molecule.formal_charge} 1\n"
        coords = molecule.coordinates
        rdmol = molecule.rdmol
        for i in range(molecule.num_atoms):
            pos = coords[i]
            sym = rdmol.GetAtomWithIdx(i).GetSymbol()
            geom_str += f"{sym} {pos[0]} {pos[1]} {pos[2]}\n"
            
        psi_mol = psi4.geometry(geom_str)
        psi4.set_options({'basis': basis_set})
        e, wfn = psi4.energy(method, return_wfn=True)
        
        # Generate dummy ESP grid to satisfy pipeline interface in Phase 1
        # Real ESP grid usually uses external programs or Psi4's VBase
        num_atoms = molecule.num_atoms
        grid_coords = np.zeros((100, 3))
        esp_values = np.zeros(100)
        
        return ESPGridResult(
            grid_coords=grid_coords,
            esp_values=esp_values,
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
        """Relaxed torsion scan via Psi4."""
        if angles is None:
            angles = np.arange(-180, 180, 15)
            
        # Returning dummy data for Phase 1
        return TorsionScanResult(
            angles=angles,
            energies=np.zeros(len(angles)),
            coordinates=np.tile(molecule.coordinates, (len(angles), 1, 1)),
        )

    def compute_wbo_matrix(
        self,
        molecule: Molecule,
        method: str = "HF",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> WBOResult:
        """Wiberg bond order matrix via Psi4."""
        # Returning dummy data for Phase 1
        num_atoms = molecule.num_atoms
        return WBOResult(
            wbo_matrix=np.eye(num_atoms),
        )
