"""
poltype.qm.psi4_backend – Psi4 QM backend adapter.

This module provides the :class:`Psi4Backend` class that wraps Psi4 for
geometry optimisation, ESP grid calculation, relaxed torsion scans, and
Wiberg–Löwdin bond-order matrix computation.

Notes on threading
------------------
The global ``psi4`` module is **not** thread-safe; all calls are serialised
through ``_PSI4_LOCK``.  Do not run two ``Psi4Backend`` instances
concurrently inside the same Python process.

Notes on CWD
------------
``psi4.oeprop`` with ``GRID_ESP`` reads ``grid.dat`` from the process CWD.
To satisfy this requirement without a permanent ``os.chdir`` the private
helper ``_in_dir`` temporarily changes CWD under the lock and restores it
on exit.
"""

from __future__ import annotations

import os
import threading
from contextlib import contextmanager
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from poltype.molecule.molecule import Molecule
from poltype.qm.backend import (
    ESPGridResult,
    OptimizationResult,
    QMBackend,
    TorsionScanResult,
    WBOResult,
)

_PSI4_LOCK = threading.Lock()

BOHR_TO_ANG: float = 0.529177210903

_VDW_RADII: dict = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98,
    "Si": 2.10, "Li": 1.82, "Na": 2.27, "K": 2.75,
}
_DEFAULT_VDW: float = 1.70


@contextmanager
def _in_dir(path: Path):
    """Context manager: cd into *path*, restore CWD on exit."""
    cwd = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd)


class Psi4Backend(QMBackend):
    """Psi4-based QM backend.

    Parameters
    ----------
    work_dir:
        Scratch directory for Psi4 input/output files.
    psi4_args:
        Reserved for future subprocess-mode use.
    num_proc:
        Number of parallel threads for Psi4.
    max_mem_gb:
        Memory limit in GB passed to Psi4.
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

    def _get_work_dir(self, molecule: Molecule) -> Path:
        wd = self.work_dir or molecule.work_dir
        wd.mkdir(parents=True, exist_ok=True)
        return wd

    def _build_geom_string(
        self,
        molecule: Molecule,
        multiplicity: int = 1,
        coords: Optional[np.ndarray] = None,
    ) -> str:
        """Return a Psi4 Cartesian geometry block (Angstroms)."""
        if coords is None:
            coords = molecule.coordinates
        rdmol = molecule.rdmol
        lines: List[str] = [f"{molecule.formal_charge} {multiplicity}"]
        for i in range(molecule.num_atoms):
            sym = rdmol.GetAtomWithIdx(i).GetSymbol()
            x, y, z = coords[i]
            lines.append(f"  {sym}  {x:.10f}  {y:.10f}  {z:.10f}")
        lines += ["  units angstrom", "  no_reorient", "  no_com"]
        return "\n".join(lines)

    def _init_psi4(self, work_dir: Path, log_name: str) -> None:
        """Reset psi4 state and point output at *work_dir/log_name*."""
        import psi4
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.set_num_threads(self.num_proc)
        psi4.set_memory(f"{self.max_mem_gb} GB")
        psi4.core.set_output_file(str(work_dir / log_name), False)
        psi4.core.IOManager.shared_object().set_default_path(str(work_dir))

    @staticmethod
    def _fibonacci_sphere(n_points: int, radius: float) -> np.ndarray:
        """Return *n_points* evenly distributed on a sphere of *radius* Å."""
        golden = (1.0 + np.sqrt(5.0)) / 2.0
        idx = np.arange(n_points, dtype=float)
        theta = np.arccos(1.0 - 2.0 * (idx + 0.5) / n_points)
        phi = 2.0 * np.pi * idx / golden
        return radius * np.column_stack([
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ])

    def _generate_esp_grid(self, molecule: Molecule) -> np.ndarray:
        """Build a multi-layer Connolly-style ESP grid around *molecule*.

        Four shells are placed at 1.4, 1.6, 1.8, and 2.0 × the Bondi VDW
        radius of each atom.  Points that fall inside *any* atom's scaled
        VDW sphere are discarded.  194 Fibonacci-lattice points per atom
        per shell are used.
        """
        coords = molecule.coordinates
        rdmol = molecule.rdmol
        n = molecule.num_atoms
        radii = np.array([
            _VDW_RADII.get(rdmol.GetAtomWithIdx(i).GetSymbol(), _DEFAULT_VDW)
            for i in range(n)
        ])
        layers: List[np.ndarray] = []
        for scale in (1.4, 1.6, 1.8, 2.0):
            scaled = radii * scale
            candidates: List[np.ndarray] = [
                self._fibonacci_sphere(194, scaled[i]) + coords[i]
                for i in range(n)
            ]
            pts = np.vstack(candidates)
            # shape: (M, N) – distance from every candidate to every atom
            dists = np.linalg.norm(
                pts[:, np.newaxis, :] - coords[np.newaxis, :, :], axis=2
            )
            outside = np.all(dists >= scaled[np.newaxis, :], axis=1)
            if outside.any():
                layers.append(pts[outside])
        return np.vstack(layers) if layers else np.empty((0, 3))

    def _get_dihedral_atoms(
        self,
        molecule: Molecule,
        rotatable_bond: Tuple[int, int],
    ) -> Tuple[int, int, int, int]:
        """Return (a, i, j, b) 0-based indices that define the dihedral."""
        i, j = rotatable_bond
        rdmol = molecule.rdmol

        def _pick_neighbor(center: int, exclude: int) -> int:
            neighbors = [
                nb.GetIdx()
                for nb in rdmol.GetAtomWithIdx(center).GetNeighbors()
                if nb.GetIdx() != exclude
            ]
            if not neighbors:
                raise ValueError(
                    f"Atom {center} has no neighbors other than {exclude}; "
                    "cannot define a dihedral."
                )
            heavy = [
                nb for nb in neighbors
                if rdmol.GetAtomWithIdx(nb).GetAtomicNum() > 1
            ]
            return heavy[0] if heavy else neighbors[0]

        return _pick_neighbor(i, j), i, j, _pick_neighbor(j, i)

    @staticmethod
    def _coords_from_psi_mol(psi_mol, n_atoms: int) -> np.ndarray:
        """Extract (n_atoms, 3) coordinates in Angstroms from a psi4 Molecule."""
        mat = psi_mol.geometry()
        return np.array(mat.np[:n_atoms, :3]) * BOHR_TO_ANG

    def optimize_geometry(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> OptimizationResult:
        """Geometry optimisation via Psi4."""
        import psi4

        work_dir = self._get_work_dir(molecule)
        multiplicity: int = kwargs.get("multiplicity", 1)

        with _PSI4_LOCK:
            self._init_psi4(work_dir, "opt.log")
            psi_mol = psi4.geometry(
                self._build_geom_string(molecule, multiplicity)
            )
            psi4.set_options({"basis": basis_set})
            try:
                e, wfn = psi4.optimize(method, return_wfn=True, molecule=psi_mol)
            except Exception as exc:
                raise RuntimeError(f"Psi4 optimization failed: {exc}") from exc

            coords = self._coords_from_psi_mol(psi_mol, molecule.num_atoms)

        return OptimizationResult(
            coordinates=coords,
            energy=float(e),
            converged=True,
            log_path=work_dir / "opt.log",
        )

    def compute_esp_grid(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "aug-cc-pVTZ",
        **kwargs,
    ) -> ESPGridResult:
        """ESP grid calculation via Psi4.

        Workflow
        --------
        1. Generate a multi-layer Connolly surface grid (Å).
        2. Run an energy calculation to obtain the wavefunction.
        3. Write ``grid.dat`` (Bohr) to *work_dir*.
        4. Call ``psi4.oeprop(wfn, 'GRID_ESP')`` – psi4 reads ``grid.dat``
           from CWD and writes ``grid_esp.dat``.  The CWD is temporarily set
           to *work_dir* for this call.
        5. Load and return ``grid_esp.dat``.
        """
        import psi4

        work_dir = self._get_work_dir(molecule)
        multiplicity: int = kwargs.get("multiplicity", 1)
        grid_ang = self._generate_esp_grid(molecule)

        with _PSI4_LOCK:
            self._init_psi4(work_dir, "esp.log")
            psi_mol = psi4.geometry(
                self._build_geom_string(molecule, multiplicity)
            )
            psi4.set_options({"basis": basis_set})
            e, wfn = psi4.energy(method, return_wfn=True, molecule=psi_mol)

            grid_bohr = grid_ang / BOHR_TO_ANG
            np.savetxt(str(work_dir / "grid.dat"), grid_bohr, fmt="%.10f")

            with _in_dir(work_dir):
                psi4.oeprop(wfn, "GRID_ESP")

            esp_values = np.loadtxt(str(work_dir / "grid_esp.dat"))

        return ESPGridResult(
            grid_coords=grid_ang,
            esp_values=esp_values,
            log_path=work_dir / "esp.log",
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

        For each requested dihedral angle a constrained geometry optimisation
        is performed using optking's ``fixed_dihedral`` keyword.  Atom indices
        in the constraint string are 1-based as required by psi4/optking.
        """
        import psi4

        if angles is None:
            angles = np.arange(-180.0, 180.0, 15.0)
        angles = np.asarray(angles, dtype=float)

        work_dir = self._get_work_dir(molecule)
        multiplicity: int = kwargs.get("multiplicity", 1)
        a1, a2, a3, a4 = self._get_dihedral_atoms(molecule, rotatable_bond)

        scan_energies: List[float] = []
        scan_coords: List[np.ndarray] = []
        log_paths: List[Path] = []

        for k, angle in enumerate(angles):
            log_name = f"torsion_{k:03d}_{angle:+.1f}deg.log"
            log_path = work_dir / log_name
            log_paths.append(log_path)

            with _PSI4_LOCK:
                self._init_psi4(work_dir, log_name)
                psi_mol = psi4.geometry(
                    self._build_geom_string(molecule, multiplicity)
                )
                psi4.set_options({
                    "basis": basis_set,
                    "optking": {
                        "fixed_dihedral": (
                            f'"{a1 + 1} {a2 + 1} {a3 + 1} {a4 + 1} {angle:.6f}"'
                        ),
                    },
                })
                try:
                    e, wfn = psi4.optimize(
                        method, return_wfn=True, molecule=psi_mol
                    )
                except Exception as exc:
                    raise RuntimeError(
                        f"Psi4 torsion scan failed at angle {angle:.1f}°: {exc}"
                    ) from exc

                coords = self._coords_from_psi_mol(psi_mol, molecule.num_atoms)

            scan_energies.append(float(e))
            scan_coords.append(coords)

        return TorsionScanResult(
            angles=angles,
            energies=np.array(scan_energies),
            coordinates=np.array(scan_coords),
            log_paths=log_paths,
        )

    def compute_wbo_matrix(
        self,
        molecule: Molecule,
        method: str = "HF",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> WBOResult:
        """Wiberg–Löwdin bond order matrix via Psi4.

        Calls ``psi4.oeprop(wfn, 'WIBERG_LOWDIN_INDICES')`` after the SCF
        and retrieves the result from ``wfn.get_array('WIBERG LOWDIN INDICES')``.
        The returned matrix has shape ``(N, N)`` where N is the number of atoms.
        """
        import psi4

        work_dir = self._get_work_dir(molecule)
        multiplicity: int = kwargs.get("multiplicity", 1)

        with _PSI4_LOCK:
            self._init_psi4(work_dir, "wbo.log")
            psi_mol = psi4.geometry(
                self._build_geom_string(molecule, multiplicity)
            )
            psi4.set_options({"basis": basis_set})
            e, wfn = psi4.energy(method, return_wfn=True, molecule=psi_mol)
            psi4.oeprop(wfn, "WIBERG_LOWDIN_INDICES")
            wbo_matrix = np.array(wfn.get_array("WIBERG LOWDIN INDICES").np)

        return WBOResult(
            wbo_matrix=wbo_matrix,
            log_path=work_dir / "wbo.log",
        )

