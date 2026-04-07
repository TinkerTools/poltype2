"""
poltype.qm.pyscf_backend – PySCF QM backend adapter.
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

# Bondi van der Waals radii in Angstrom, used for ESP grid generation.
_VDW_RADII: dict[str, float] = {
    "H": 1.20, "He": 1.40,
    "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52,
    "F": 1.47, "Ne": 1.54,
    "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80,
    "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Br": 1.85, "I": 1.98,
}
_VDW_DEFAULT = 1.80


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

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _build_pyscf_mol(self, molecule: Molecule, basis_set: str):
        """Construct and build a ``pyscf.gto.Mole`` from *molecule*."""
        from pyscf import gto, lib

        lib.num_threads(self.num_proc)
        rdmol = molecule.rdmol
        conf = rdmol.GetConformer()
        atom_str = "; ".join(
            f"{rdmol.GetAtomWithIdx(i).GetSymbol()} "
            f"{conf.GetAtomPosition(i).x:.10f} "
            f"{conf.GetAtomPosition(i).y:.10f} "
            f"{conf.GetAtomPosition(i).z:.10f}"
            for i in range(molecule.num_atoms)
        )
        mol = gto.Mole()
        mol.atom = atom_str
        mol.basis = basis_set
        mol.charge = molecule.formal_charge
        mol.spin = 0
        mol.max_memory = self.max_mem_gb * 1024
        mol.verbose = 0
        mol.build()
        return mol

    def _get_mf_and_dm(self, pyscf_mol, method: str):
        """Run the requested QM method and return ``(mf, dm_ao)``."""
        from pyscf import scf, mp, dft

        method_upper = method.upper()

        if method_upper in ("HF", "RHF"):
            mf = scf.RHF(pyscf_mol)
            mf.kernel()
            dm = mf.make_rdm1()
            return mf, dm

        if method_upper == "MP2":
            mf = scf.RHF(pyscf_mol)
            mf.kernel()
            pt = mp.MP2(mf)
            pt.kernel()
            dm_mo = pt.make_rdm1()
            # Transform 1-RDM from MO to AO basis.
            mo = mf.mo_coeff
            dm = mo @ dm_mo @ mo.T
            return mf, dm

        # Treat everything else as a DFT functional.
        mf = dft.RKS(pyscf_mol)
        mf.xc = method
        mf.kernel()
        dm = mf.make_rdm1()
        return mf, dm

    def _build_mf_for_opt(self, pyscf_mol, method: str):
        """Return an MF/correlated object suitable for geomeTRIC optimisation."""
        from pyscf import scf, mp, dft

        method_upper = method.upper()
        if method_upper in ("HF", "RHF"):
            return scf.RHF(pyscf_mol)
        if method_upper == "MP2":
            mf = scf.RHF(pyscf_mol)
            return mp.MP2(mf)
        mf = dft.RKS(pyscf_mol)
        mf.xc = method
        return mf

    def _build_esp_grid(
        self,
        molecule: Molecule,
        pyscf_mol,
        spacing: float = 0.3,
        outer_margin: float = 2.8,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(grid_ang, grid_bohr)`` for a CHELPG-style ESP grid.

        A uniform cubic lattice covers the molecular bounding box plus
        *outer_margin* Å on each side.  Points inside any atomic VDW sphere
        are discarded.
        """
        from pyscf import lib

        coords_ang = molecule.coordinates  # (N, 3), Angstrom
        rdmol = molecule.rdmol
        vdw = np.array(
            [
                _VDW_RADII.get(rdmol.GetAtomWithIdx(i).GetSymbol(), _VDW_DEFAULT)
                for i in range(molecule.num_atoms)
            ]
        )

        lo = coords_ang.min(axis=0) - outer_margin
        hi = coords_ang.max(axis=0) + outer_margin
        x = np.arange(lo[0], hi[0] + spacing, spacing)
        y = np.arange(lo[1], hi[1] + spacing, spacing)
        z = np.arange(lo[2], hi[2] + spacing, spacing)
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        grid = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

        # diff[p, a, :] = grid_point_p - atom_a; shape (npts, natom, 3)
        diff = grid[:, None, :] - coords_ang[None, :, :]
        dists = np.linalg.norm(diff, axis=2)  # (npts, natom)

        # Keep points outside every atomic VDW sphere.
        outside_vdw = (dists - vdw[None, :]).min(axis=1) > 0.0
        grid_ang = grid[outside_vdw]
        grid_bohr = grid_ang / lib.param.BOHR
        return grid_ang, grid_bohr

    def _find_dihedral_atoms(
        self, molecule: Molecule, rotatable_bond: Tuple[int, int]
    ) -> Tuple[int, int, int, int]:
        """Return ``(i, j, k, l)`` 0-based indices for the dihedral spanning *rotatable_bond*."""
        j, k = rotatable_bond
        rdmol = molecule.rdmol
        i_candidates = [
            nb.GetIdx()
            for nb in rdmol.GetAtomWithIdx(j).GetNeighbors()
            if nb.GetIdx() != k
        ]
        l_candidates = [
            nb.GetIdx()
            for nb in rdmol.GetAtomWithIdx(k).GetNeighbors()
            if nb.GetIdx() != j
        ]
        if not i_candidates or not l_candidates:
            raise ValueError(
                f"Cannot determine dihedral atoms for bond ({j}, {k}): "
                f"i_candidates={i_candidates}, l_candidates={l_candidates}"
            )
        return i_candidates[0], j, k, l_candidates[0]

    # ------------------------------------------------------------------
    # QMBackend interface
    # ------------------------------------------------------------------

    def optimize_geometry(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> OptimizationResult:
        """Geometry optimisation via PySCF + geomeTRIC."""
        try:
            from pyscf import lib
            from pyscf.geomopt.geometric_solver import optimize
        except ImportError as e:
            raise ImportError(
                "PySCF and geomeTRIC are required for optimize_geometry. "
                "Install them with `pip install pyscf geometric`."
            ) from e

        pyscf_mol = self._build_pyscf_mol(molecule, basis_set)
        mf = self._build_mf_for_opt(pyscf_mol, method)
        mf_opt = optimize(mf)

        final_mol = mf_opt
        if not hasattr(final_mol, "atom_coords"):
            final_mol = getattr(mf_opt, "mol", mf_opt)
        
        if not hasattr(final_mol, "atom_coords"):
            raise RuntimeError(
                f"Could not extract optimised coordinates from {type(mf_opt)}"
            )

        final_coords = final_mol.atom_coords() * lib.param.BOHR
        
        # Extract energy. Some correlated methods in PySCF don't have e_tot 
        # on the optimizer object but keep it in the underlying mf/correlated object.
        final_energy = 0.0
        try:
            final_energy = float(getattr(final_mol, "e_tot", 0.0))
        except (TypeError, ValueError):
            pass

        if final_energy == 0.0:
            try:
                if hasattr(mf_opt, "e_tot"):
                    final_energy = float(mf_opt.e_tot)
                elif hasattr(mf, "e_tot"):
                    final_energy = float(mf.e_tot)
            except (TypeError, ValueError):
                pass

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
        """Compute the molecular electrostatic potential on a CHELPG-style grid.

        The ESP at each grid point r is::

            V(r) = sum_A Z_A / |r - R_A|
                   - sum_{mu,nu} dm_{mu nu} * (mu nu | 1/r | s_r)

        where s_r is a unit-charge s-type Gaussian at r, evaluated via
        ``pyscf.df.incore.aux_e2`` with the ``int3c2e`` kernel.
        Units: atomic units (Hartree/e).
        """
        try:
            from pyscf import gto, df, lib
        except ImportError as e:
            raise ImportError("PySCF is required for compute_esp_grid.") from e

        pyscf_mol = self._build_pyscf_mol(molecule, basis_set)
        _mf, dm = self._get_mf_and_dm(pyscf_mol, method)

        grid_ang, grid_bohr = self._build_esp_grid(molecule, pyscf_mol)
        npts = len(grid_ang)

        # Nuclear electrostatic contribution: sum_A Z_A / |r - R_A|
        v_nuc = np.zeros(npts)
        for idx in range(pyscf_mol.natm):
            r_atom = pyscf_mol.atom_coord(idx)  # Bohr
            z_atom = pyscf_mol.atom_charge(idx)
            dr = grid_bohr - r_atom
            v_nuc += z_atom / np.linalg.norm(dr, axis=1)

        # Electronic contribution via 3-centre 2-electron integrals.
        # eri3c[mu, nu, p] = (mu nu | 1/r_1p | s_p)
        fakemol = gto.fakemol_for_charges(grid_bohr)
        eri3c = df.incore.aux_e2(pyscf_mol, fakemol, intor="int3c2e")
        v_elec = -np.einsum("ijp,ij->p", eri3c, dm)

        return ESPGridResult(
            grid_coords=grid_ang,
            esp_values=v_nuc + v_elec,
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
        """Relaxed torsion energy scan via PySCF + geomeTRIC constrained optimisation.

        For each target dihedral angle a constrained geometry optimisation is
        performed.  The geometry from step n is used as the starting point for
        step n+1 to improve convergence along the scan.  Failed steps store
        ``nan`` energies and fall back to the previous geometry.
        """
        try:
            from pyscf import lib
            from pyscf.geomopt.geometric_solver import optimize
        except ImportError as e:
            raise ImportError(
                "PySCF and geomeTRIC are required for torsion_scan."
            ) from e

        if angles is None:
            angles = np.arange(-180, 180, 15, dtype=float)
        angles = np.asarray(angles, dtype=float)

        i, j, k, l = self._find_dihedral_atoms(molecule, rotatable_bond)
        # geomeTRIC constraint files use 1-based atom indices.
        i1, j1, k1, l1 = i + 1, j + 1, k + 1, l + 1

        all_energies: list[float] = []
        all_coords: list[np.ndarray] = []
        current_molecule = molecule

        for angle_deg in angles:
            constraint_str = (
                f"$set\n"
                f"dihedral {i1} {j1} {k1} {l1} {angle_deg:.4f}\n"
                f"$end\n"
            )
            pyscf_mol = self._build_pyscf_mol(current_molecule, basis_set)
            mf = self._build_mf_for_opt(pyscf_mol, method)
            try:
                mf_opt = optimize(mf, constraints=constraint_str)
                energy = float(getattr(mf_opt, "e_tot", float("nan")))
                opt_mol = getattr(mf_opt, "mol", None)
                if opt_mol is not None and hasattr(opt_mol, "atom_coords"):
                    coords_ang = opt_mol.atom_coords() * lib.param.BOHR
                else:
                    coords_ang = current_molecule.coordinates.copy()
                current_molecule = molecule.with_updated_coordinates(coords_ang)
            except Exception:
                energy = float("nan")
                coords_ang = current_molecule.coordinates.copy()

            all_energies.append(energy)
            all_coords.append(coords_ang)

        return TorsionScanResult(
            angles=angles,
            energies=np.array(all_energies),
            coordinates=np.array(all_coords),
        )

    def compute_wbo_matrix(
        self,
        molecule: Molecule,
        method: str = "HF",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> WBOResult:
        """Compute the Mayer bond-order matrix via PySCF.

        The Mayer bond order between atoms A and B is::

            BO(A, B) = sum_{mu in A, nu in B} (D S)_{mu nu} (D S)_{nu mu}

        where D is the AO density matrix and S is the AO overlap matrix.
        """
        pyscf_mol = self._build_pyscf_mol(molecule, basis_set)
        _mf, dm = self._get_mf_and_dm(pyscf_mol, method)

        S = pyscf_mol.intor_symmetric("int1e_ovlp")
        DS = dm @ S

        natm = pyscf_mol.natm
        ao_labels = pyscf_mol.ao_labels(fmt=False)
        ao_atom = np.array([lbl[0] for lbl in ao_labels])

        wbo = np.zeros((natm, natm))
        for A in range(natm):
            idx_A = np.where(ao_atom == A)[0]
            for B in range(A + 1, natm):
                idx_B = np.where(ao_atom == B)[0]
                block_AB = DS[np.ix_(idx_A, idx_B)]
                block_BA = DS[np.ix_(idx_B, idx_A)]
                bo = float(np.sum(block_AB * block_BA.T))
                wbo[A, B] = bo
                wbo[B, A] = bo

        return WBOResult(wbo_matrix=wbo)

