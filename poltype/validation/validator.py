"""
poltype.validation.validator – abstract and concrete parameter validators.

Every validator implements the :class:`Validator` ABC.  The concrete
:class:`EnergyValidator` compares torsion-scan QM energies against MM
estimates to flag poor fits.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Dict, List, Optional

import numpy as np

from poltype.validation.validation_result import (
    ValidationRecord,
    ValidationReport,
    ValidationSeverity,
)

if TYPE_CHECKING:
    from poltype.pipeline.context import PipelineContext


# ---------------------------------------------------------------------------
# Abstract base
# ---------------------------------------------------------------------------


class Validator(ABC):
    """Abstract base class for parameter validators.

    Subclasses must implement :meth:`validate` which inspects the
    pipeline context (artifacts, molecule, config) and returns a
    :class:`ValidationReport`.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable validator identifier."""

    @abstractmethod
    def validate(self, context: "PipelineContext") -> ValidationReport:
        """Run validation checks and return a report.

        Parameters
        ----------
        context:
            The pipeline context containing artifacts from upstream
            stages.

        Returns
        -------
        ValidationReport
        """

    def __repr__(self) -> str:
        return f"{type(self).__name__}()"


# ---------------------------------------------------------------------------
# Concrete: energy-based validation
# ---------------------------------------------------------------------------


class EnergyValidator(Validator):
    """Compare torsion-scan QM energies against MM-estimated energies.

    The validator examines the ``torsion_scans`` artifact and checks:

    1. **Energy convergence** – the optimisation energy is finite and
       the geometry converged.
    2. **Torsion energy range** – the span of torsion energies is
       within a configurable tolerance.
    3. **Torsion RMSD** – when MM torsion estimates are available, the
       root-mean-square deviation between QM and MM energies is below
       a threshold.

    Parameters
    ----------
    energy_range_threshold:
        Maximum allowed range (max − min) in the torsion energy
        profile, in kcal/mol.  Default ``50.0``.
    rmsd_threshold:
        Maximum allowed RMSD between QM and MM torsion profiles,
        in kcal/mol.  Default ``3.0``.
    """

    HARTREE_TO_KCALMOL = 627.509474

    def __init__(
        self,
        energy_range_threshold: float = 50.0,
        rmsd_threshold: float = 3.0,
    ) -> None:
        self._energy_range_threshold = energy_range_threshold
        self._rmsd_threshold = rmsd_threshold

    @property
    def name(self) -> str:
        return "energy"

    @property
    def energy_range_threshold(self) -> float:
        """Maximum allowed torsion-energy range (kcal/mol)."""
        return self._energy_range_threshold

    @property
    def rmsd_threshold(self) -> float:
        """Maximum allowed QM-vs-MM RMSD (kcal/mol)."""
        return self._rmsd_threshold

    def validate(self, context: "PipelineContext") -> ValidationReport:
        """Run energy validation checks.

        Returns
        -------
        ValidationReport
            A report with one record per check.
        """
        records: List[ValidationRecord] = []

        mol_name = context.molecule.name if context.molecule else "unknown"

        # --- check 1: geometry optimisation energy ---
        opt_result = context.get_artifact("opt_result")
        records.append(self._check_opt_energy(opt_result))

        # --- check 2: torsion energy range ---
        torsion_scans = context.get_artifact("torsion_scans")
        records.extend(self._check_torsion_range(torsion_scans))

        # --- check 3: QM vs MM RMSD ---
        mm_torsion_energies = context.get_artifact("mm_torsion_energies")
        records.extend(
            self._check_torsion_rmsd(torsion_scans, mm_torsion_energies)
        )

        return ValidationReport(
            records=records,
            molecule_name=mol_name,
            validation_method=self.name,
            metadata={
                "energy_range_threshold": self._energy_range_threshold,
                "rmsd_threshold": self._rmsd_threshold,
            },
        )

    # -- internal check helpers -----------------------------------------------

    def _check_opt_energy(
        self, opt_result: Optional[Any],
    ) -> ValidationRecord:
        """Validate geometry optimisation energy."""
        if opt_result is None:
            return ValidationRecord(
                check_name="energy_convergence",
                passed=True,
                severity=ValidationSeverity.INFO,
                message="No optimisation result; skipping energy convergence check",
            )

        converged = getattr(opt_result, "converged", True)
        energy = getattr(opt_result, "energy", None)

        if energy is not None and not np.isfinite(energy):
            return ValidationRecord(
                check_name="energy_convergence",
                passed=False,
                severity=ValidationSeverity.ERROR,
                message=f"Optimisation energy is not finite: {energy}",
                details={"energy": float(energy), "converged": converged},
            )

        if not converged:
            return ValidationRecord(
                check_name="energy_convergence",
                passed=False,
                severity=ValidationSeverity.WARNING,
                message="Geometry optimisation did not converge",
                details={"energy": float(energy) if energy is not None else None, "converged": False},
            )

        return ValidationRecord(
            check_name="energy_convergence",
            passed=True,
            severity=ValidationSeverity.INFO,
            message="Geometry optimisation converged",
            details={"energy": float(energy) if energy is not None else None, "converged": True},
        )

    def _check_torsion_range(
        self, torsion_scans: Optional[Dict],
    ) -> List[ValidationRecord]:
        """Validate that torsion energy ranges are within tolerance."""
        if not torsion_scans:
            return [
                ValidationRecord(
                    check_name="torsion_energy_range",
                    passed=True,
                    severity=ValidationSeverity.INFO,
                    message="No torsion scans to validate",
                ),
            ]

        records: List[ValidationRecord] = []
        for bond, scan in torsion_scans.items():
            energies_hartree = getattr(scan, "energies", None)
            if energies_hartree is None or len(energies_hartree) == 0:
                records.append(
                    ValidationRecord(
                        check_name="torsion_energy_range",
                        passed=True,
                        severity=ValidationSeverity.INFO,
                        message=f"Bond {bond}: no energy data",
                        details={"bond": bond},
                    )
                )
                continue

            energies_kcal = (
                np.asarray(energies_hartree, dtype=float) * self.HARTREE_TO_KCALMOL
            )
            energy_range = float(np.ptp(energies_kcal))
            passed = energy_range <= self._energy_range_threshold

            records.append(
                ValidationRecord(
                    check_name="torsion_energy_range",
                    passed=passed,
                    severity=ValidationSeverity.WARNING if not passed else ValidationSeverity.INFO,
                    message=(
                        f"Bond {bond}: energy range {energy_range:.2f} kcal/mol"
                        + ("" if passed else f" exceeds threshold {self._energy_range_threshold}")
                    ),
                    details={
                        "bond": bond,
                        "energy_range_kcal": energy_range,
                        "threshold": self._energy_range_threshold,
                    },
                )
            )

        return records

    def _check_torsion_rmsd(
        self,
        torsion_scans: Optional[Dict],
        mm_energies: Optional[Dict],
    ) -> List[ValidationRecord]:
        """Compare QM and MM torsion profiles when both are available."""
        if not torsion_scans or not mm_energies:
            return [
                ValidationRecord(
                    check_name="torsion_rmsd",
                    passed=True,
                    severity=ValidationSeverity.INFO,
                    message="QM/MM torsion comparison not available",
                ),
            ]

        records: List[ValidationRecord] = []
        for bond, scan in torsion_scans.items():
            mm_scan = mm_energies.get(bond)
            if mm_scan is None:
                continue

            qm = np.asarray(getattr(scan, "energies", []), dtype=float)
            mm = np.asarray(mm_scan, dtype=float)

            if qm.size == 0 or mm.size == 0 or qm.shape != mm.shape:
                continue

            # Convert to kcal/mol relative energies
            qm_rel = (qm - qm.min()) * self.HARTREE_TO_KCALMOL
            mm_rel = mm - mm.min()
            rmsd = float(np.sqrt(np.mean((qm_rel - mm_rel) ** 2)))
            passed = rmsd <= self._rmsd_threshold

            records.append(
                ValidationRecord(
                    check_name="torsion_rmsd",
                    passed=passed,
                    severity=ValidationSeverity.WARNING if not passed else ValidationSeverity.INFO,
                    message=(
                        f"Bond {bond}: QM/MM RMSD {rmsd:.3f} kcal/mol"
                        + ("" if passed else f" exceeds threshold {self._rmsd_threshold}")
                    ),
                    details={
                        "bond": bond,
                        "rmsd_kcal": rmsd,
                        "threshold": self._rmsd_threshold,
                    },
                )
            )

        return records if records else [
            ValidationRecord(
                check_name="torsion_rmsd",
                passed=True,
                severity=ValidationSeverity.INFO,
                message="No matching QM/MM torsion data for comparison",
            ),
        ]
