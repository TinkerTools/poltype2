"""
poltype.validation.validation_result – data classes for validation outcomes.

A :class:`ValidationRecord` holds the outcome of a single validation
check, and a :class:`ValidationReport` aggregates all records produced
during the validation stage.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Any, Dict, List


class ValidationSeverity(enum.Enum):
    """Severity of a validation finding."""

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


@dataclass(frozen=True)
class ValidationRecord:
    """A single validation finding.

    Attributes
    ----------
    check_name:
        Short identifier for the validation check (e.g.
        ``"energy_convergence"``).
    passed:
        Whether the check passed.
    severity:
        How serious a failure is — ``INFO``, ``WARNING``, or ``ERROR``.
    message:
        Human-readable description of the finding.
    details:
        Optional key/value map with numeric or diagnostic details.
    """

    check_name: str
    passed: bool
    severity: ValidationSeverity = ValidationSeverity.INFO
    message: str = ""
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ValidationReport:
    """Aggregated report produced by the validation stage.

    Attributes
    ----------
    records:
        Individual validation findings.
    molecule_name:
        Name of the molecule being validated.
    validation_method:
        Label describing the validator(s) used.
    metadata:
        Extra information (e.g. thresholds, method names).
    """

    records: List[ValidationRecord] = field(default_factory=list)
    molecule_name: str = ""
    validation_method: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    # -- convenience properties --

    @property
    def num_records(self) -> int:
        """Total number of validation records."""
        return len(self.records)

    @property
    def num_passed(self) -> int:
        """Number of checks that passed."""
        return sum(1 for r in self.records if r.passed)

    @property
    def num_failed(self) -> int:
        """Number of checks that did **not** pass."""
        return sum(1 for r in self.records if not r.passed)

    @property
    def all_passed(self) -> bool:
        """``True`` when every check passed (or there are no checks)."""
        return all(r.passed for r in self.records)

    @property
    def has_errors(self) -> bool:
        """``True`` when at least one record has ``ERROR`` severity and did not pass."""
        return any(
            not r.passed and r.severity is ValidationSeverity.ERROR
            for r in self.records
        )

    @property
    def warnings(self) -> List[ValidationRecord]:
        """Failed records with ``WARNING`` severity."""
        return [
            r for r in self.records
            if not r.passed and r.severity is ValidationSeverity.WARNING
        ]

    @property
    def errors(self) -> List[ValidationRecord]:
        """Failed records with ``ERROR`` severity."""
        return [
            r for r in self.records
            if not r.passed and r.severity is ValidationSeverity.ERROR
        ]
