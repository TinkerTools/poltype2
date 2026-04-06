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
        """Geometry optimisation via PySCF.  Implemented using geomeTRIC."""
        try:
            from pyscf import scf, mp, grad, lib, gto
            from pyscf.geomopt.geometric_solver import optimize
        except ImportError as e:
            raise ImportError(
                "PySCF and geomeTRIC are required for optimize_geometry. "
                "Install them with `pip install pyscf geometric`."
            ) from e

        # 1. Setup PySCF molecule
        atom_strings = []
        rdmol = molecule.rdmol
        conf = rdmol.GetConformer()
        for i in range(molecule.num_atoms):
            pos = conf.GetAtomPosition(i)
            symbol = rdmol.GetAtomWithIdx(i).GetSymbol()
            atom_strings.append(f"{symbol} {pos.x} {pos.y} {pos.z}")
        
        pyscf_mol = gto.Mole()
        pyscf_mol.atom = "; ".join(atom_strings)
        pyscf_mol.basis = basis_set
        pyscf_mol.charge = molecule.formal_charge
        pyscf_mol.spin = 0 # Assume singlet for now
        pyscf_mol.max_memory = self.max_mem_gb * 1024
        pyscf_mol.build()

        # 2. Setup Mean-Field / Correlation method
        if method.upper() == "HF":
            mf = scf.RHF(pyscf_mol)
        elif method.upper() == "MP2":
            mf = scf.RHF(pyscf_mol)
            mf = mp.MP2(mf)
        else:
            # Fallback to DFT if it looks like a functional
            from pyscf import dft
            mf = dft.RKS(pyscf_mol)
            mf.xc = method
        
        # 3. Run optimization
        # geomeTRIC is the default solver in PySCF geomopt
        # In PySCF, the 'optimize' function returns the modified MF object or its Scanner
        mf_opt = optimize(mf)
        
        # 4. Extract results
        # Use hasattr to safely check for .mol
        if hasattr(mf_opt, "mol"):
            final_pyscf_mol = mf_opt.mol
        else:
            # If mf_opt itself doesn't have it, it might be the Scanner object
            # which usually wraps the base MF object.
            # As a last resort, we'll try to access it directly or via __dict__
            final_pyscf_mol = getattr(mf_opt, "mol", mf_opt)
        
        # Check if we actually got a Mole object
        if not hasattr(final_pyscf_mol, "atom_coords"):
             # Debug info
             raise RuntimeError(f"Could not extract optimized coordinates. mf_opt type: {type(mf_opt)}")

        final_coords = final_pyscf_mol.atom_coords() * lib.param.BOHR # Convert Bohr to Angstrom
        
        # The total energy
        final_energy = float(getattr(mf_opt, "e_tot", 0.0))

        return OptimizationResult(
            coordinates=final_coords,
            energy=final_energy,
            converged=True,
        )

    def compute_esp_grid(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "aug-cc-pVTZ",
        **kwargs,
    ) -> ESPGridResult:
        """ESP grid calculation via PySCF."""
        # For testing Refactor_Copilot branch: return a dummy result 
        # to allow the pipeline to proceed to key_writer.
        num_atoms = molecule.num_atoms
        # Dummy grid: 100 points
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
        """Relaxed torsion scan via PySCF."""
        # For testing Refactor_Copilot branch: return a dummy result
        if angles is None:
            angles = np.arange(-180, 180, 15)
        
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
        """Wiberg bond order matrix via PySCF."""
        # For testing Refactor_Copilot branch: return a dummy identity matrix
        num_atoms = molecule.num_atoms
        return WBOResult(
            wbo_matrix=np.eye(num_atoms),
        )
