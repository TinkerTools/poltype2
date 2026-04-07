"""
poltype.qm.gaussian_backend – Gaussian QM backend adapter.

Implements geometry optimisation, ESP grid calculation, relaxed torsion
scan, and Wiberg bond order matrix calculation by writing Gaussian input
files, running the ``g16`` (or configurable) executable, and parsing the
resulting log files.
"""

from __future__ import annotations

import copy
import re
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
from rdkit.Chem import AllChem

from poltype.molecule.molecule import Molecule
from poltype.qm.backend import (
    DMAResult,
    ESPGridResult,
    OptimizationResult,
    QMBackend,
    TorsionScanResult,
    WBOResult,
)

_BOHR_TO_ANG = 0.529177210903


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

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _resolve_work_dir(self, molecule: Molecule) -> Path:
        wd = self.work_dir if self.work_dir is not None else molecule.work_dir
        wd.mkdir(parents=True, exist_ok=True)
        return wd

    def _link0(self) -> str:
        return f"%NProcShared={self.num_proc}\n%Mem={self.max_mem_gb}GB\n"

    def _coord_block(self, molecule: Molecule) -> str:
        """Return ``charge multiplicity`` line followed by Cartesian coordinates."""
        rdmol = molecule.rdmol
        conf = rdmol.GetConformer(0)
        lines = [f"{molecule.formal_charge} 1"]
        for i in range(molecule.num_atoms):
            pos = conf.GetAtomPosition(i)
            sym = rdmol.GetAtomWithIdx(i).GetSymbol()
            lines.append(
                f" {sym:<2}  {pos.x:>14.8f}  {pos.y:>14.8f}  {pos.z:>14.8f}"
            )
        return "\n".join(lines)

    def _run_gaussian(self, inp_path: Path) -> Path:
        """Execute Gaussian on *inp_path* (cwd = parent dir) and return the log path."""
        log_path = inp_path.with_suffix(".log")
        result = subprocess.run(
            [self.gaus_exe, inp_path.name],
            cwd=inp_path.parent,
            capture_output=True,
            text=True,
        )
        if not log_path.exists():
            raise RuntimeError(
                f"Gaussian did not produce a log file for {inp_path}.\n"
                f"stdout: {result.stdout[:1000]}\n"
                f"stderr: {result.stderr[:1000]}"
            )
        return log_path

    # ------------------------------------------------------------------
    # Parsing helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_final_energy(log_text: str) -> float:
        """Extract the last reported energy from a Gaussian log.

        Tries MP2 (EUMP2), then DFT/HF SCF Done, in that order.
        """
        mp2 = re.findall(r"EUMP2\s*=\s*([-\d.DE+]+)", log_text)
        if mp2:
            return float(mp2[-1].replace("D", "E"))
        scf = re.findall(r"SCF Done:\s*E\([^)]+\)\s*=\s*([-\d.]+)", log_text)
        if scf:
            return float(scf[-1])
        raise ValueError("Could not parse final energy from Gaussian log.")

    @staticmethod
    def _parse_orientation_block(log_text: str) -> np.ndarray:
        """Parse the last Standard orientation (or Input orientation) block.

        Returns an ``(N, 3)`` float array of Cartesian coordinates in Angstroms.
        """
        for label in ("Standard orientation", "Input orientation"):
            blocks = re.findall(
                rf"{label}:.*?-{{5,}}\n(.*?)\n\s*-{{5,}}",
                log_text,
                re.DOTALL,
            )
            if blocks:
                coords: List[List[float]] = []
                for line in blocks[-1].strip().splitlines():
                    parts = line.split()
                    if len(parts) >= 6:
                        try:
                            coords.append(
                                [float(parts[3]), float(parts[4]), float(parts[5])]
                            )
                        except ValueError:
                            continue
                if coords:
                    return np.array(coords)
        raise ValueError("No orientation block found in Gaussian log.")

    @staticmethod
    def _parse_esp_grid(log_text: str) -> Tuple[np.ndarray, np.ndarray]:
        """Extract ESP grid coordinates (Angstroms) and potentials (a.u.) from the log.

        Gaussian writes the raw grid when IOp(6/33=2) is active; coordinates
        are in Bohr in the log and are converted here.
        """
        section = re.search(
            r"Electrostatic Properties Using The .*?Density.*?"
            r"-{5,}\n(.*?)\n\s*-{5,}",
            log_text,
            re.DOTALL,
        )
        if section is None:
            section = re.search(
                r"Grid points\s+X\s+Y\s+Z\s+Potential.*?\n(.*?)\n\s*-{5,}",
                log_text,
                re.DOTALL,
            )
        if section is None:
            raise ValueError(
                "Cannot locate ESP grid section in Gaussian log. "
                "Ensure the route includes Pop=MK IOp(6/33=2)."
            )
        raw = section.group(1)
        grid_coords: List[List[float]] = []
        esp_values: List[float] = []
        for line in raw.strip().splitlines():
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                x, y, z, v = (
                    float(parts[0]),
                    float(parts[1]),
                    float(parts[2]),
                    float(parts[3]),
                )
            except ValueError:
                continue
            grid_coords.append(
                [x * _BOHR_TO_ANG, y * _BOHR_TO_ANG, z * _BOHR_TO_ANG]
            )
            esp_values.append(v)
        return np.array(grid_coords), np.array(esp_values)

    @staticmethod
    def _parse_wbo_matrix(log_text: str, num_atoms: int) -> np.ndarray:
        """Parse the Wiberg bond index matrix from Gaussian NBO output.

        The matrix is printed in column-blocked form (one block per page);
        this routine accumulates all blocks into a single (N, N) array.
        """
        wbo_match = re.search(
            r"Wiberg bond index matrix in the NAO basis:(.*?)(?:\n\s*\n|\Z)",
            log_text,
            re.DOTALL,
        )
        if wbo_match is None:
            raise ValueError(
                "Wiberg bond index matrix not found in Gaussian log. "
                "Ensure the route includes Pop=NBORead and $NBO BNDIDX $END."
            )
        raw_lines = [
            ln for ln in wbo_match.group(1).strip().splitlines() if ln.strip()
        ]
        matrix = np.zeros((num_atoms, num_atoms))
        col_indices: List[int] = []
        for line in raw_lines:
            parts = line.split()
            # Column-header lines: every token parses as an integer
            try:
                col_indices = [int(p) - 1 for p in parts]
                continue
            except ValueError:
                pass
            # Data line: first token is 1-based row index, rest are floats
            if not col_indices:
                continue
            try:
                row_idx = int(parts[0]) - 1
                values = [float(v) for v in parts[1:]]
            except (ValueError, IndexError):
                continue
            for j, v in zip(col_indices, values):
                if 0 <= row_idx < num_atoms and 0 <= j < num_atoms:
                    matrix[row_idx, j] = v
        return matrix

    @staticmethod
    def _dihedral_atoms(
        molecule: Molecule, rotatable_bond: Tuple[int, int]
    ) -> Tuple[int, int, int, int]:
        """Return four 0-based atom indices forming the scanned dihedral."""
        rdmol = molecule.rdmol
        i_bond, j_bond = rotatable_bond
        nbrs_i = [
            n.GetIdx()
            for n in rdmol.GetAtomWithIdx(i_bond).GetNeighbors()
            if n.GetIdx() != j_bond
        ]
        nbrs_j = [
            n.GetIdx()
            for n in rdmol.GetAtomWithIdx(j_bond).GetNeighbors()
            if n.GetIdx() != i_bond
        ]
        if not nbrs_i or not nbrs_j:
            raise ValueError(
                f"Cannot construct a dihedral for bond {rotatable_bond}: "
                "one atom has no neighbours other than the bond partner."
            )
        return nbrs_i[0], i_bond, j_bond, nbrs_j[0]

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
        """Geometry optimisation via Gaussian.

        Writes a ``<name>_opt.gjf`` input, runs Gaussian, then parses the
        final standard-orientation coordinates and energy from the log.

        Keyword arguments
        -----------------
        max_cycles : int
            Maximum optimisation cycles (default 100).
        convergence : str
            Gaussian convergence keyword, e.g. ``"Tight"`` (default: none).
        """
        wd = self._resolve_work_dir(molecule)
        stem = molecule.name or "mol"
        inp_path = wd / f"{stem}_opt.gjf"

        max_cycles = int(kwargs.get("max_cycles", 100))
        convergence = kwargs.get("convergence", "")
        conv_str = f",{convergence}" if convergence else ""
        route = (
            f"#P {method}/{basis_set} "
            f"Opt=(MaxCycles={max_cycles}{conv_str}) NoSymm"
        )
        content = (
            f"{self._link0()}"
            f"%Chk={stem}_opt.chk\n"
            f"{route}\n\n"
            f"{stem} geometry optimisation\n\n"
            f"{self._coord_block(molecule)}\n\n"
        )
        inp_path.write_text(content)
        log_path = self._run_gaussian(inp_path)
        log_text = log_path.read_text()

        converged = bool(
            re.search(
                r"Optimization completed\.|Stationary point found\.", log_text
            )
        )
        try:
            energy = self._parse_final_energy(log_text)
        except ValueError as exc:
            raise RuntimeError(
                f"Energy parse failed for {log_path}: {exc}"
            ) from exc

        coords = self._parse_orientation_block(log_text)
        return OptimizationResult(
            coordinates=coords,
            energy=energy,
            converged=converged,
            log_path=log_path,
        )

    def compute_dma(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "6-311G**",
        **kwargs,
    ) -> DMAResult:
        """DMA single-point via Gaussian to produce an fchk file.

        Writes a Gaussian input that runs a single-point calculation and
        produces a checkpoint file.  The checkpoint is then converted to a
        formatted checkpoint (fchk) using ``formchk``.  The resulting fchk
        file is consumed by GDMA.

        Keyword arguments
        -----------------
        formchk_exe : str
            Path to the ``formchk`` executable (default ``"formchk"``).
        """
        wd = self._resolve_work_dir(molecule)
        stem = molecule.name or "mol"
        inp_path = wd / f"{stem}_dma.gjf"
        chk_name = f"{stem}_dma.chk"
        chk_path = wd / chk_name

        density_kw = " Density=MP2" if method.upper() == "MP2" else ""
        route = (
            f"#P {method}/{basis_set}{density_kw} "
            f"SP NoSymm"
        )
        content = (
            f"{self._link0()}"
            f"%Chk={chk_name}\n"
            f"{route}\n\n"
            f"{stem} DMA single-point\n\n"
            f"{self._coord_block(molecule)}\n\n"
        )
        inp_path.write_text(content)
        log_path = self._run_gaussian(inp_path)
        log_text = log_path.read_text()

        try:
            energy = self._parse_final_energy(log_text)
        except ValueError as exc:
            raise RuntimeError(
                f"Energy parse failed for DMA calculation {log_path}: {exc}"
            ) from exc

        # Convert chk → fchk
        fchk_path = wd / f"{stem}_dma.fchk"
        formchk_exe = kwargs.get("formchk_exe", "formchk")
        result = subprocess.run(
            [formchk_exe, str(chk_path), str(fchk_path)],
            cwd=str(wd),
            capture_output=True,
            text=True,
        )
        if not fchk_path.exists():
            raise RuntimeError(
                f"formchk did not produce {fchk_path}.\n"
                f"stdout: {result.stdout[:1000]}\n"
                f"stderr: {result.stderr[:1000]}"
            )

        return DMAResult(
            fchk_path=fchk_path,
            energy=energy,
            log_path=log_path,
        )

    def compute_esp_grid(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "aug-cc-pVTZ",
        **kwargs,
    ) -> ESPGridResult:
        """ESP single-point with Merz-Kollman grid output via Gaussian.

        Uses Pop=MK IOp(6/33=2,6/41=10,6/42=17) to write raw grid coordinates
        and potentials to the log file.  For MP2, Density=MP2 is added so the
        ESP reflects the correlated density rather than the HF reference.
        """
        wd = self._resolve_work_dir(molecule)
        stem = molecule.name or "mol"
        inp_path = wd / f"{stem}_esp.gjf"

        density_kw = " Density=MP2" if method.upper() == "MP2" else ""
        route = (
            f"#P {method}/{basis_set}{density_kw} "
            f"Pop=MK IOp(6/33=2,6/41=10,6/42=17) NoSymm"
        )
        content = (
            f"{self._link0()}"
            f"%Chk={stem}_esp.chk\n"
            f"{route}\n\n"
            f"{stem} ESP single-point\n\n"
            f"{self._coord_block(molecule)}\n\n"
        )
        inp_path.write_text(content)
        log_path = self._run_gaussian(inp_path)
        log_text = log_path.read_text()

        grid_coords, esp_values = self._parse_esp_grid(log_text)
        return ESPGridResult(
            grid_coords=grid_coords,
            esp_values=esp_values,
            log_path=log_path,
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
        """Relaxed torsion scan via Gaussian (one constrained optimisation per angle).

        For each target dihedral angle the starting geometry is pre-set using
        RDKit, then a Gaussian Opt=(ModRedundant) job is launched with the
        dihedral frozen (``F``).  Running one job per angle supports arbitrary
        non-uniform grids and avoids corner cases in Gaussian's built-in scan.

        Keyword arguments
        -----------------
        max_cycles : int
            Maximum optimisation cycles per scan point (default 50).
        """
        if angles is None:
            angles = np.arange(-180, 180, 15, dtype=float)

        wd = self._resolve_work_dir(molecule)
        stem = molecule.name or "mol"
        max_cycles = int(kwargs.get("max_cycles", 50))

        a1_0, a2_0, a3_0, a4_0 = self._dihedral_atoms(molecule, rotatable_bond)
        # Gaussian uses 1-based atom indices in the ModRedundant section
        g_a1, g_a2, g_a3, g_a4 = a1_0 + 1, a2_0 + 1, a3_0 + 1, a4_0 + 1

        scan_energies: List[float] = []
        scan_coords: List[np.ndarray] = []
        log_paths: List[Path] = []

        for k, angle in enumerate(angles):
            # Pre-set the dihedral in a deep copy so the F constraint is exact
            rdmol_k = copy.deepcopy(molecule.rdmol)
            AllChem.SetDihedralDeg(
                rdmol_k.GetConformer(0),
                a1_0, a2_0, a3_0, a4_0,
                float(angle),
            )
            mol_k = Molecule(rdmol_k, work_dir=molecule.work_dir, name=molecule.name)

            inp_path = wd / f"{stem}_scan_{k:03d}.gjf"
            route = (
                f"#P {method}/{basis_set} "
                f"Opt=(ModRedundant,MaxCycles={max_cycles}) NoSymm"
            )
            frozen_line = f"D {g_a1} {g_a2} {g_a3} {g_a4} F"
            content = (
                f"{self._link0()}"
                f"%Chk={stem}_scan_{k:03d}.chk\n"
                f"{route}\n\n"
                f"{stem} torsion scan step {k} angle {angle:.2f}\n\n"
                f"{self._coord_block(mol_k)}\n\n"
                f"{frozen_line}\n\n"
            )
            inp_path.write_text(content)
            log_path = self._run_gaussian(inp_path)
            log_text = log_path.read_text()
            log_paths.append(log_path)

            try:
                energy = self._parse_final_energy(log_text)
            except ValueError as exc:
                raise RuntimeError(
                    f"Energy parse failed at angle {angle:.2f} ({log_path}): {exc}"
                ) from exc

            coords = self._parse_orientation_block(log_text)
            scan_energies.append(energy)
            scan_coords.append(coords)

        return TorsionScanResult(
            angles=np.array(angles, dtype=float),
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
        """Single-point NBO analysis to obtain the Wiberg bond order matrix.

        The route uses Pop=NBORead; the NBO input section requests
        $NBO BNDIDX $END which causes Gaussian/NBO to print the full Wiberg
        bond index matrix in the NAO basis to the log file.
        """
        wd = self._resolve_work_dir(molecule)
        stem = molecule.name or "mol"
        inp_path = wd / f"{stem}_wbo.gjf"

        route = f"#P {method}/{basis_set} Pop=NBORead NoSymm"
        content = (
            f"{self._link0()}"
            f"%Chk={stem}_wbo.chk\n"
            f"{route}\n\n"
            f"{stem} WBO single-point\n\n"
            f"{self._coord_block(molecule)}\n\n"
            "$NBO BNDIDX $END\n\n"
        )
        inp_path.write_text(content)
        log_path = self._run_gaussian(inp_path)
        log_text = log_path.read_text()

        wbo_matrix = self._parse_wbo_matrix(log_text, molecule.num_atoms)
        return WBOResult(
            wbo_matrix=wbo_matrix,
            log_path=log_path,
        )

