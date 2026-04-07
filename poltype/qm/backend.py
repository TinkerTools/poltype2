"""
poltype.qm.backend – abstract QM backend interface.

All QM backends (Psi4, Gaussian, PySCF, ANI, …) must implement the
:class:`QMBackend` abstract base class defined here.  Pipeline stages
interact **only** with this interface and never inspect which concrete
backend is in use.

Result types are plain :mod:`dataclasses` so they are serialisable and
easy to inspect in tests.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

from poltype.molecule.molecule import Molecule


# ---------------------------------------------------------------------------
# Typed result containers
# ---------------------------------------------------------------------------


@dataclass
class OptimizationResult:
    """Outcome of a QM geometry optimisation.

    Attributes
    ----------
    coordinates:
        Optimised Cartesian coordinates, shape ``(N, 3)``, Angstroms.
    energy:
        Final SCF / MP2 / DFT energy in Hartree.
    converged:
        Whether the optimisation converged to the requested threshold.
    log_path:
        Path to the QM log / output file (may be ``None`` if not written).
    """

    coordinates: np.ndarray
    energy: float
    converged: bool = True
    log_path: Optional[Path] = None


@dataclass
class ESPGridResult:
    """Outcome of an electrostatic potential grid calculation.

    Attributes
    ----------
    grid_coords:
        Cartesian coordinates of the ESP grid points, shape ``(M, 3)``,
        Angstroms.
    esp_values:
        Electrostatic potential at each grid point, shape ``(M,)``,
        in atomic units (Hartree / e).
    log_path:
        Path to the QM log / output file.
    """

    grid_coords: np.ndarray
    esp_values: np.ndarray
    log_path: Optional[Path] = None


@dataclass
class TorsionScanResult:
    """Outcome of a relaxed torsion energy scan.

    Attributes
    ----------
    angles:
        Dihedral angles scanned, in degrees, shape ``(K,)``.
    energies:
        QM energies at each angle, shape ``(K,)``, in Hartree.
    coordinates:
        Optimised Cartesian coordinates at each angle, shape ``(K, N, 3)``.
    log_paths:
        Per-point QM log file paths.
    """

    angles: np.ndarray
    energies: np.ndarray
    coordinates: np.ndarray = field(
        default_factory=lambda: np.empty((0, 0, 3))
    )
    log_paths: list = field(default_factory=list)


@dataclass
class DMAResult:
    """Outcome of a DMA (Distributed Multipole Analysis) QM calculation.

    This calculation is typically run at MP2/6-311G** to produce a
    wavefunction that GDMA can analyse for atomic multipoles.

    Attributes
    ----------
    fchk_path:
        Path to the formatted checkpoint file produced by the calculation.
        This file is consumed by GDMA.
    energy:
        Final QM energy in Hartree.
    log_path:
        Path to the QM log / output file.
    """

    fchk_path: Path
    energy: float = 0.0
    log_path: Optional[Path] = None


@dataclass
class WBOResult:
    """Outcome of a Wiberg Bond Order matrix calculation.

    Attributes
    ----------
    wbo_matrix:
        Wiberg bond order matrix, shape ``(N, N)``.
    log_path:
        Path to the QM log / output file.
    """

    wbo_matrix: np.ndarray
    log_path: Optional[Path] = None


# ---------------------------------------------------------------------------
# Abstract backend
# ---------------------------------------------------------------------------


class QMBackend(ABC):
    """Abstract base class for QM calculation backends.

    Concrete subclasses implement the four core QM operations required by
    the Poltype2 pipeline.  All methods accept a :class:`~poltype.molecule.molecule.Molecule`
    as input and return a typed result object.

    All file I/O should be performed relative to ``molecule.work_dir``
    (absolute path); backends must **never** call ``os.chdir``.

    Parameters
    ----------
    work_dir:
        Base directory for scratch / output files.  If ``None``, each
        method should use ``molecule.work_dir``.
    """

    def __init__(self, work_dir: Optional[Path] = None) -> None:
        self.work_dir = work_dir

    # ------------------------------------------------------------------
    # Abstract methods – every backend must implement these four
    # ------------------------------------------------------------------

    @abstractmethod
    def optimize_geometry(
        self,
        molecule: Molecule,
        method: str,
        basis_set: str,
        **kwargs,
    ) -> OptimizationResult:
        """Run a QM geometry optimisation.

        Parameters
        ----------
        molecule:
            Input molecule with a starting 3-D conformer.
        method:
            QM method string (e.g. ``"MP2"``, ``"wB97X-D"``).
        basis_set:
            Basis set string (e.g. ``"6-31G*"``).
        **kwargs:
            Backend-specific keyword arguments (e.g. ``max_cycles``,
            ``convergence``).

        Returns
        -------
        OptimizationResult
        """

    @abstractmethod
    def compute_esp_grid(
        self,
        molecule: Molecule,
        method: str,
        basis_set: str,
        **kwargs,
    ) -> ESPGridResult:
        """Compute the electrostatic potential on a grid around *molecule*.

        Parameters
        ----------
        molecule:
            Input molecule at the target geometry.
        method:
            QM method string.
        basis_set:
            Basis set string.
        **kwargs:
            Backend-specific keyword arguments.

        Returns
        -------
        ESPGridResult
        """

    @abstractmethod
    def torsion_scan(
        self,
        molecule: Molecule,
        rotatable_bond: Tuple[int, int],
        method: str,
        basis_set: str,
        angles: Optional[np.ndarray] = None,
        **kwargs,
    ) -> TorsionScanResult:
        """Perform a relaxed torsion energy scan.

        Parameters
        ----------
        molecule:
            Input molecule.
        rotatable_bond:
            ``(i, j)`` 0-based atom indices defining the scanned dihedral.
        method:
            QM method string.
        basis_set:
            Basis set string.
        angles:
            Array of dihedral angles to scan (degrees).  If ``None``, a
            default grid from −180 to 175 in 15° steps is used.
        **kwargs:
            Backend-specific keyword arguments.

        Returns
        -------
        TorsionScanResult
        """

    @abstractmethod
    def compute_dma(
        self,
        molecule: Molecule,
        method: str,
        basis_set: str,
        **kwargs,
    ) -> DMAResult:
        """Run a QM single-point to produce a formatted checkpoint (fchk).

        The resulting fchk file is consumed by the GDMA program to
        extract distributed multipoles.  Typically called with
        ``method="MP2"`` and ``basis_set="6-311G**"``.

        Parameters
        ----------
        molecule:
            Input molecule at the optimised geometry.
        method:
            QM method string (e.g. ``"MP2"``).
        basis_set:
            Basis set string (e.g. ``"6-311G**"``).
        **kwargs:
            Backend-specific keyword arguments.

        Returns
        -------
        DMAResult
        """

    @abstractmethod
    def compute_wbo_matrix(
        self,
        molecule: Molecule,
        method: str,
        basis_set: str,
        **kwargs,
    ) -> WBOResult:
        """Compute the Wiberg Bond Order (WBO) matrix.

        Parameters
        ----------
        molecule:
            Input molecule.
        method:
            QM method string.
        basis_set:
            Basis set string.
        **kwargs:
            Backend-specific keyword arguments.

        Returns
        -------
        WBOResult
        """

    # ------------------------------------------------------------------
    # Optional / convenience methods
    # ------------------------------------------------------------------

    @property
    def name(self) -> str:
        """Human-readable name of this backend (e.g. ``"psi4"``)."""
        return type(self).__name__

    def __repr__(self) -> str:
        return f"{self.name}(work_dir={self.work_dir})"
