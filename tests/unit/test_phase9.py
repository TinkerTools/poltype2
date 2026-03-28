"""
tests/unit/test_phase9.py – unit tests for Phase 9: validation & quality assurance.

Tests cover:
- ValidationSeverity enum values.
- ValidationRecord dataclass (frozen, attributes, severity).
- ValidationReport dataclass (properties, aggregation, filters).
- Validator ABC cannot be instantiated.
- EnergyValidator: opt energy check, torsion range check, RMSD check,
  threshold configuration, missing artifacts, edge cases.
- ValidationStage: execute, skip logic, no molecule, default validator,
  injected validator, artifact storage.
- Updated pipeline factory: 10 stages, stage order.
- Full pipeline integration: validation stage completes.
- Dry-run mode: validation stage skips.
- Finalization collects validation_report.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pytest

from poltype.config.schema import PoltypeConfig, QMConfig
from poltype.molecule.molecule import Molecule
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.factory import build_default_pipeline
from poltype.pipeline.stage import StageResult, StageStatus
from poltype.pipeline.stages.validation import ValidationStage
from poltype.qm.backend import (
    ESPGridResult,
    OptimizationResult,
    QMBackend,
    TorsionScanResult,
    WBOResult,
)
from poltype.validation.validation_result import (
    ValidationRecord,
    ValidationReport,
    ValidationSeverity,
)
from poltype.validation.validator import EnergyValidator, Validator


# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------


class MockBackend(QMBackend):
    """In-memory mock backend for integration tests."""

    @property
    def name(self) -> str:
        return "mock"

    def optimize_geometry(self, molecule, method="MP2", basis_set="6-31G*", **kw):
        return OptimizationResult(
            coordinates=molecule.coordinates,
            energy=-76.123456,
            converged=True,
        )

    def compute_esp_grid(self, molecule, method="MP2", basis_set="aug-cc-pVTZ", **kw):
        n_points = 50
        return ESPGridResult(
            grid_coords=np.random.randn(n_points, 3),
            esp_values=np.random.randn(n_points),
        )

    def torsion_scan(
        self, molecule, rotatable_bond=(0, 1), method="wB97X-D",
        basis_set="6-31G*", angles=None, **kw,
    ):
        scan_angles = angles if angles is not None else np.arange(-180, 180, 15.0)
        return TorsionScanResult(
            angles=scan_angles,
            energies=np.zeros_like(scan_angles),
        )

    def compute_wbo_matrix(self, molecule, method="HF", basis_set="6-31G*", **kw):
        n = molecule.num_atoms
        return WBOResult(wbo_matrix=np.eye(n))


def _make_ethanol() -> Molecule:
    """Ethanol – has C, O, H and a rotatable bond for richer tests."""
    return Molecule.from_smiles("CCO", name="ethanol")


def _make_methane() -> Molecule:
    """Methane – simplest hydrocarbon for quick tests."""
    return Molecule.from_smiles("C", name="methane")


def _make_context(
    tmp_path: Path,
    molecule: Optional[Molecule] = None,
    backend: Optional[QMBackend] = None,
    dry_run: bool = False,
    skip_validation: bool = False,
) -> PipelineContext:
    config = PoltypeConfig(
        dry_run=dry_run,
        skip_validation=skip_validation,
    )
    mol = molecule or _make_ethanol()
    return PipelineContext(
        config=config,
        molecule=mol,
        backend=backend or MockBackend(),
        work_dir=tmp_path,
    )


# ===================================================================
# A. ValidationSeverity tests
# ===================================================================


class TestValidationSeverity:
    def test_info_value(self):
        assert ValidationSeverity.INFO.value == "info"

    def test_warning_value(self):
        assert ValidationSeverity.WARNING.value == "warning"

    def test_error_value(self):
        assert ValidationSeverity.ERROR.value == "error"

    def test_all_members(self):
        assert set(ValidationSeverity) == {
            ValidationSeverity.INFO,
            ValidationSeverity.WARNING,
            ValidationSeverity.ERROR,
        }


# ===================================================================
# B. ValidationRecord tests
# ===================================================================


class TestValidationRecord:
    def test_creation_minimal(self):
        rec = ValidationRecord(check_name="test_check", passed=True)
        assert rec.check_name == "test_check"
        assert rec.passed is True
        assert rec.severity is ValidationSeverity.INFO
        assert rec.message == ""
        assert rec.details == {}

    def test_creation_full(self):
        rec = ValidationRecord(
            check_name="energy_convergence",
            passed=False,
            severity=ValidationSeverity.ERROR,
            message="Energy is not finite",
            details={"energy": float("inf")},
        )
        assert rec.check_name == "energy_convergence"
        assert rec.passed is False
        assert rec.severity is ValidationSeverity.ERROR
        assert rec.message == "Energy is not finite"
        assert "energy" in rec.details

    def test_frozen(self):
        rec = ValidationRecord(check_name="test", passed=True)
        with pytest.raises(AttributeError):
            rec.passed = False  # type: ignore[misc]

    def test_equality(self):
        a = ValidationRecord(check_name="x", passed=True)
        b = ValidationRecord(check_name="x", passed=True)
        assert a == b

    def test_inequality(self):
        a = ValidationRecord(check_name="x", passed=True)
        b = ValidationRecord(check_name="x", passed=False)
        assert a != b


# ===================================================================
# C. ValidationReport tests
# ===================================================================


class TestValidationReport:
    def test_empty_report(self):
        report = ValidationReport()
        assert report.num_records == 0
        assert report.num_passed == 0
        assert report.num_failed == 0
        assert report.all_passed is True
        assert report.has_errors is False
        assert report.warnings == []
        assert report.errors == []

    def test_all_passed(self):
        report = ValidationReport(
            records=[
                ValidationRecord(check_name="a", passed=True),
                ValidationRecord(check_name="b", passed=True),
            ],
        )
        assert report.num_records == 2
        assert report.num_passed == 2
        assert report.num_failed == 0
        assert report.all_passed is True

    def test_mixed_results(self):
        report = ValidationReport(
            records=[
                ValidationRecord(check_name="a", passed=True),
                ValidationRecord(
                    check_name="b", passed=False,
                    severity=ValidationSeverity.WARNING,
                ),
                ValidationRecord(
                    check_name="c", passed=False,
                    severity=ValidationSeverity.ERROR,
                ),
            ],
        )
        assert report.num_records == 3
        assert report.num_passed == 1
        assert report.num_failed == 2
        assert report.all_passed is False
        assert report.has_errors is True
        assert len(report.warnings) == 1
        assert len(report.errors) == 1

    def test_warnings_only(self):
        report = ValidationReport(
            records=[
                ValidationRecord(
                    check_name="a", passed=False,
                    severity=ValidationSeverity.WARNING,
                ),
            ],
        )
        assert report.has_errors is False
        assert len(report.warnings) == 1

    def test_molecule_name(self):
        report = ValidationReport(molecule_name="ethanol")
        assert report.molecule_name == "ethanol"

    def test_validation_method(self):
        report = ValidationReport(validation_method="energy")
        assert report.validation_method == "energy"

    def test_metadata(self):
        report = ValidationReport(metadata={"threshold": 3.0})
        assert report.metadata["threshold"] == 3.0

    def test_passed_filter_excludes_passed_records(self):
        """warnings and errors properties only return failed records."""
        rec = ValidationRecord(
            check_name="a", passed=True,
            severity=ValidationSeverity.WARNING,
        )
        report = ValidationReport(records=[rec])
        assert report.warnings == []


# ===================================================================
# D. Validator ABC tests
# ===================================================================


class TestValidatorABC:
    def test_cannot_instantiate(self):
        with pytest.raises(TypeError):
            Validator()  # type: ignore[abstract]

    def test_concrete_must_implement_validate(self):
        class PartialValidator(Validator):
            @property
            def name(self):
                return "partial"

        with pytest.raises(TypeError):
            PartialValidator()  # type: ignore[abstract]

    def test_concrete_must_implement_name(self):
        class PartialValidator(Validator):
            def validate(self, context):
                return ValidationReport()

        with pytest.raises(TypeError):
            PartialValidator()  # type: ignore[abstract]


# ===================================================================
# E. EnergyValidator tests
# ===================================================================


class TestEnergyValidator:
    def test_name(self):
        v = EnergyValidator()
        assert v.name == "energy"

    def test_repr(self):
        v = EnergyValidator()
        assert "EnergyValidator" in repr(v)

    def test_default_thresholds(self):
        v = EnergyValidator()
        assert v.energy_range_threshold == 50.0
        assert v.rmsd_threshold == 3.0

    def test_custom_thresholds(self):
        v = EnergyValidator(energy_range_threshold=100.0, rmsd_threshold=5.0)
        assert v.energy_range_threshold == 100.0
        assert v.rmsd_threshold == 5.0


class TestEnergyValidatorOptCheck:
    def test_no_opt_result(self, tmp_path):
        """When no opt_result artifact, the check passes (info)."""
        v = EnergyValidator()
        ctx = _make_context(tmp_path)
        report = v.validate(ctx)
        convergence = [r for r in report.records if r.check_name == "energy_convergence"]
        assert len(convergence) == 1
        assert convergence[0].passed is True

    def test_converged_opt(self, tmp_path):
        """A converged optimisation passes the energy_convergence check."""
        v = EnergyValidator()
        ctx = _make_context(tmp_path)
        ctx.set_artifact("opt_result", OptimizationResult(
            coordinates=np.zeros((3, 3)), energy=-76.0, converged=True,
        ))
        report = v.validate(ctx)
        convergence = [r for r in report.records if r.check_name == "energy_convergence"]
        assert convergence[0].passed is True
        assert convergence[0].details["converged"] is True

    def test_unconverged_opt(self, tmp_path):
        """An unconverged optimisation triggers a warning."""
        v = EnergyValidator()
        ctx = _make_context(tmp_path)
        ctx.set_artifact("opt_result", OptimizationResult(
            coordinates=np.zeros((3, 3)), energy=-76.0, converged=False,
        ))
        report = v.validate(ctx)
        convergence = [r for r in report.records if r.check_name == "energy_convergence"]
        assert convergence[0].passed is False
        assert convergence[0].severity is ValidationSeverity.WARNING

    def test_nan_energy_fails(self, tmp_path):
        """A NaN optimisation energy triggers an error."""
        v = EnergyValidator()
        ctx = _make_context(tmp_path)
        ctx.set_artifact("opt_result", OptimizationResult(
            coordinates=np.zeros((3, 3)), energy=float("nan"), converged=True,
        ))
        report = v.validate(ctx)
        convergence = [r for r in report.records if r.check_name == "energy_convergence"]
        assert convergence[0].passed is False
        assert convergence[0].severity is ValidationSeverity.ERROR

    def test_inf_energy_fails(self, tmp_path):
        """An infinite optimisation energy triggers an error."""
        v = EnergyValidator()
        ctx = _make_context(tmp_path)
        ctx.set_artifact("opt_result", OptimizationResult(
            coordinates=np.zeros((3, 3)), energy=float("inf"), converged=True,
        ))
        report = v.validate(ctx)
        convergence = [r for r in report.records if r.check_name == "energy_convergence"]
        assert convergence[0].passed is False
        assert convergence[0].severity is ValidationSeverity.ERROR


class TestEnergyValidatorTorsionRange:
    def test_no_torsion_scans(self, tmp_path):
        """When no torsion_scans artifact, the check passes (info)."""
        v = EnergyValidator()
        ctx = _make_context(tmp_path)
        report = v.validate(ctx)
        range_checks = [r for r in report.records if r.check_name == "torsion_energy_range"]
        assert len(range_checks) == 1
        assert range_checks[0].passed is True

    def test_flat_energy_profile_passes(self, tmp_path):
        """A flat torsion energy profile is within any positive threshold."""
        v = EnergyValidator(energy_range_threshold=50.0)
        ctx = _make_context(tmp_path)
        angles = np.arange(-180, 180, 15.0)
        scan = TorsionScanResult(
            angles=angles,
            energies=np.full_like(angles, -76.0),
        )
        ctx.set_artifact("torsion_scans", {(0, 1): scan})
        report = v.validate(ctx)
        range_checks = [r for r in report.records if r.check_name == "torsion_energy_range"]
        assert range_checks[0].passed is True

    def test_large_range_fails(self, tmp_path):
        """A large energy range exceeding the threshold triggers a warning."""
        v = EnergyValidator(energy_range_threshold=10.0)
        ctx = _make_context(tmp_path)
        angles = np.arange(-180, 180, 15.0)
        # 0.1 Hartree range ≈ 62.75 kcal/mol
        energies = np.linspace(-76.0, -75.9, len(angles))
        scan = TorsionScanResult(angles=angles, energies=energies)
        ctx.set_artifact("torsion_scans", {(0, 1): scan})
        report = v.validate(ctx)
        range_checks = [r for r in report.records if r.check_name == "torsion_energy_range"]
        assert range_checks[0].passed is False
        assert range_checks[0].severity is ValidationSeverity.WARNING

    def test_multiple_bonds(self, tmp_path):
        """Each bond gets its own torsion_energy_range record."""
        v = EnergyValidator()
        ctx = _make_context(tmp_path)
        angles = np.arange(-180, 180, 15.0)
        scan_a = TorsionScanResult(angles=angles, energies=np.zeros_like(angles))
        scan_b = TorsionScanResult(angles=angles, energies=np.zeros_like(angles))
        ctx.set_artifact("torsion_scans", {(0, 1): scan_a, (2, 3): scan_b})
        report = v.validate(ctx)
        range_checks = [r for r in report.records if r.check_name == "torsion_energy_range"]
        assert len(range_checks) == 2


class TestEnergyValidatorRMSD:
    def test_no_mm_energies(self, tmp_path):
        """Without MM energies, the RMSD check passes as info."""
        v = EnergyValidator()
        ctx = _make_context(tmp_path)
        report = v.validate(ctx)
        rmsd_checks = [r for r in report.records if r.check_name == "torsion_rmsd"]
        assert len(rmsd_checks) == 1
        assert rmsd_checks[0].passed is True

    def test_perfect_match(self, tmp_path):
        """When QM and MM perfectly match, RMSD is 0 and check passes."""
        v = EnergyValidator(rmsd_threshold=1.0)
        ctx = _make_context(tmp_path)
        angles = np.arange(-180, 180, 15.0)
        energies = np.zeros_like(angles)
        scan = TorsionScanResult(angles=angles, energies=energies)
        ctx.set_artifact("torsion_scans", {(0, 1): scan})
        # MM energies in kcal/mol (relative) — all zeros match QM
        ctx.set_artifact("mm_torsion_energies", {(0, 1): np.zeros_like(angles)})
        report = v.validate(ctx)
        rmsd_checks = [r for r in report.records if r.check_name == "torsion_rmsd"]
        assert rmsd_checks[0].passed is True
        assert rmsd_checks[0].details["rmsd_kcal"] == pytest.approx(0.0, abs=1e-10)

    def test_large_rmsd_fails(self, tmp_path):
        """When RMSD exceeds threshold, the check fails."""
        v = EnergyValidator(rmsd_threshold=0.1)
        ctx = _make_context(tmp_path)
        angles = np.arange(-180, 180, 15.0)
        # QM: 0.01 Hartree range ≈ 6.275 kcal/mol
        qm_energies = np.sin(np.radians(angles)) * 0.01
        scan = TorsionScanResult(angles=angles, energies=qm_energies)
        ctx.set_artifact("torsion_scans", {(0, 1): scan})
        # MM: flat (all zeros kcal/mol)
        ctx.set_artifact("mm_torsion_energies", {(0, 1): np.zeros_like(angles)})
        report = v.validate(ctx)
        rmsd_checks = [r for r in report.records if r.check_name == "torsion_rmsd"]
        assert rmsd_checks[0].passed is False

    def test_mismatched_bonds_skipped(self, tmp_path):
        """When MM energies have different bond keys, no RMSD check is generated."""
        v = EnergyValidator()
        ctx = _make_context(tmp_path)
        angles = np.arange(-180, 180, 15.0)
        scan = TorsionScanResult(angles=angles, energies=np.zeros_like(angles))
        ctx.set_artifact("torsion_scans", {(0, 1): scan})
        ctx.set_artifact("mm_torsion_energies", {(2, 3): np.zeros_like(angles)})
        report = v.validate(ctx)
        rmsd_checks = [r for r in report.records if r.check_name == "torsion_rmsd"]
        assert all(r.passed for r in rmsd_checks)


class TestEnergyValidatorReport:
    def test_report_molecule_name(self, tmp_path):
        v = EnergyValidator()
        ctx = _make_context(tmp_path)
        report = v.validate(ctx)
        assert report.molecule_name == "ethanol"

    def test_report_method(self, tmp_path):
        v = EnergyValidator()
        ctx = _make_context(tmp_path)
        report = v.validate(ctx)
        assert report.validation_method == "energy"

    def test_report_metadata_contains_thresholds(self, tmp_path):
        v = EnergyValidator(energy_range_threshold=42.0, rmsd_threshold=2.5)
        ctx = _make_context(tmp_path)
        report = v.validate(ctx)
        assert report.metadata["energy_range_threshold"] == 42.0
        assert report.metadata["rmsd_threshold"] == 2.5

    def test_no_molecule_uses_unknown(self, tmp_path):
        v = EnergyValidator()
        config = PoltypeConfig()
        ctx = PipelineContext(config=config, work_dir=tmp_path)
        report = v.validate(ctx)
        assert report.molecule_name == "unknown"


# ===================================================================
# F. ValidationStage tests
# ===================================================================


class TestValidationStage:
    def test_name(self):
        stage = ValidationStage()
        assert stage.name == "validation"

    def test_execute_no_molecule_fails(self, tmp_path):
        stage = ValidationStage()
        config = PoltypeConfig()
        ctx = PipelineContext(config=config, work_dir=tmp_path)
        result = stage.execute(ctx)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_execute_with_molecule_completes(self, tmp_path):
        stage = ValidationStage()
        ctx = _make_context(tmp_path)
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED

    def test_execute_produces_validation_report(self, tmp_path):
        stage = ValidationStage()
        ctx = _make_context(tmp_path)
        stage.execute(ctx)
        report = ctx.get_artifact("validation_report")
        assert report is not None
        assert isinstance(report, ValidationReport)

    def test_execute_default_validator(self, tmp_path):
        """When no validators injected, EnergyValidator is used."""
        stage = ValidationStage()
        ctx = _make_context(tmp_path)
        stage.execute(ctx)
        report = ctx.get_artifact("validation_report")
        assert "energy" in report.validation_method

    def test_injected_validators(self, tmp_path):
        """Custom validators can be injected."""

        class AlwaysPassValidator(Validator):
            @property
            def name(self):
                return "always_pass"

            def validate(self, context):
                return ValidationReport(
                    records=[
                        ValidationRecord(check_name="dummy", passed=True),
                    ],
                )

        stage = ValidationStage(validators=[AlwaysPassValidator()])
        ctx = _make_context(tmp_path)
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        report = ctx.get_artifact("validation_report")
        assert report.validation_method == "always_pass"
        assert report.all_passed is True

    def test_validators_property(self):
        stage = ValidationStage()
        assert stage.validators == []

    def test_validators_property_with_injected(self):
        v = EnergyValidator()
        stage = ValidationStage(validators=[v])
        assert len(stage.validators) == 1

    def test_result_message_format(self, tmp_path):
        stage = ValidationStage()
        ctx = _make_context(tmp_path)
        result = stage.execute(ctx)
        assert "passed" in result.message
        assert "failed" in result.message

    def test_artifacts_contain_validation_report(self, tmp_path):
        stage = ValidationStage()
        ctx = _make_context(tmp_path)
        result = stage.execute(ctx)
        assert "validation_report" in result.artifacts


class TestValidationStageSkip:
    def test_should_skip_when_dry_run(self, tmp_path):
        stage = ValidationStage()
        ctx = _make_context(tmp_path, dry_run=True)
        assert stage.should_skip(ctx) is True

    def test_should_skip_when_skip_validation(self, tmp_path):
        stage = ValidationStage()
        ctx = _make_context(tmp_path, skip_validation=True)
        assert stage.should_skip(ctx) is True

    def test_should_not_skip_by_default(self, tmp_path):
        stage = ValidationStage()
        ctx = _make_context(tmp_path)
        assert stage.should_skip(ctx) is False


# ===================================================================
# G. Updated pipeline factory tests
# ===================================================================


class TestPipelineFactoryWithValidation:
    def test_ten_stages(self):
        runner = build_default_pipeline()
        assert len(runner._stages) == 10

    def test_stage_names_in_order(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert names == [
            "input_preparation",
            "geometry_optimization",
            "esp_fitting",
            "multipole",
            "atom_typing",
            "database_match",
            "fragmentation",
            "torsion_fitting",
            "validation",
            "finalization",
        ]

    def test_validation_stage_present(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert "validation" in names

    def test_validation_after_torsion(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert names.index("torsion_fitting") < names.index("validation")

    def test_validation_before_finalization(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert names.index("validation") < names.index("finalization")


# ===================================================================
# H. Full pipeline integration tests
# ===================================================================


class TestFullPipelineWithValidation:
    def test_validation_runs_in_full_pipeline(self, tmp_path):
        """ValidationStage completes in a full pipeline run."""
        mol = _make_ethanol()
        config = PoltypeConfig(
            small_molecule_smarts_to_tinker_class="/nonexistent_smarts.txt"
        )
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)
        assert "validation" in ctx.stage_results
        assert ctx.stage_results["validation"].status is StageStatus.COMPLETED

    def test_validation_report_available(self, tmp_path):
        """validation_report artifact is available after validation runs."""
        mol = _make_ethanol()
        config = PoltypeConfig(
            small_molecule_smarts_to_tinker_class="/nonexistent_smarts.txt"
        )
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)
        result = ctx.get_artifact("validation_report")
        assert result is not None
        assert isinstance(result, ValidationReport)

    def test_all_ten_stages_complete(self, tmp_path):
        """All 10 stages complete in a full pipeline run."""
        mol = _make_ethanol()
        config = PoltypeConfig(
            small_molecule_smarts_to_tinker_class="/nonexistent_smarts.txt"
        )
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)
        assert len(ctx.stage_results) == 10
        for sr in ctx.stage_results.values():
            assert sr.status is StageStatus.COMPLETED

    def test_finalization_collects_validation_report(self, tmp_path):
        """Finalization stage includes validation in summary."""
        mol = _make_ethanol()
        config = PoltypeConfig(
            small_molecule_smarts_to_tinker_class="/nonexistent_smarts.txt"
        )
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)
        summary = ctx.get_artifact("final_summary")
        assert summary is not None
        assert "validation_report" in summary["collected_artifacts"]


class TestDryRunWithValidation:
    def test_validation_skips_in_dry_run(self, tmp_path):
        stage = ValidationStage()
        ctx = _make_context(tmp_path, dry_run=True)
        assert stage.should_skip(ctx) is True

    def test_dry_run_full_pipeline_skips_validation(self, tmp_path):
        mol = _make_ethanol()
        config = PoltypeConfig(dry_run=True)
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)
        assert ctx.stage_results["validation"].status is StageStatus.SKIPPED
        # input_prep and finalization still run
        assert ctx.stage_results["input_preparation"].status is StageStatus.COMPLETED
        assert ctx.stage_results["finalization"].status is StageStatus.COMPLETED

    def test_skip_validation_flag_skips(self, tmp_path):
        mol = _make_ethanol()
        config = PoltypeConfig(skip_validation=True)
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)
        assert ctx.stage_results["validation"].status is StageStatus.SKIPPED


# ===================================================================
# I. Multiple validators
# ===================================================================


class TestMultipleValidators:
    def test_two_validators(self, tmp_path):
        """Two validators produce combined records."""

        class PassValidator(Validator):
            @property
            def name(self):
                return "pass_validator"

            def validate(self, context):
                return ValidationReport(
                    records=[ValidationRecord(check_name="pass_check", passed=True)],
                )

        class FailValidator(Validator):
            @property
            def name(self):
                return "fail_validator"

            def validate(self, context):
                return ValidationReport(
                    records=[
                        ValidationRecord(
                            check_name="fail_check", passed=False,
                            severity=ValidationSeverity.WARNING,
                        ),
                    ],
                )

        stage = ValidationStage(validators=[PassValidator(), FailValidator()])
        ctx = _make_context(tmp_path)
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        report = ctx.get_artifact("validation_report")
        assert report.num_records == 2
        assert report.num_passed == 1
        assert report.num_failed == 1
        assert "pass_validator, fail_validator" in report.validation_method


# ===================================================================
# J. Import tests
# ===================================================================


class TestModuleImports:
    def test_import_validation_module(self):
        from poltype.validation import (
            EnergyValidator,
            ValidationRecord,
            ValidationReport,
            ValidationSeverity,
            Validator,
        )
        assert Validator is not None
        assert EnergyValidator is not None
        assert ValidationRecord is not None
        assert ValidationReport is not None
        assert ValidationSeverity is not None

    def test_import_stage_from_stages_init(self):
        from poltype.pipeline.stages import ValidationStage
        assert ValidationStage is not None

    def test_import_stage_directly(self):
        from poltype.pipeline.stages.validation import ValidationStage
        assert ValidationStage is not None


# ===================================================================
# K. Config field tests
# ===================================================================


class TestConfigSkipValidation:
    def test_default_false(self):
        config = PoltypeConfig()
        assert config.skip_validation is False

    def test_set_true(self):
        config = PoltypeConfig(skip_validation=True)
        assert config.skip_validation is True
