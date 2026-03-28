"""
tests/unit/test_phase4.py – unit tests for Phase 4: stage implementations.

Tests cover:
- GeometryOptimizationStage: delegates to backend, updates molecule,
  handles convergence failure and backend errors.
- ESPFittingStage: delegates to backend, stores ESP artifact.
- MultipoleStage: consumes ESP artifact, produces multipole frames.
- TorsionFittingStage: iterates rotatable bonds, collects scan results.
- FinalizationStage: assembles summary from upstream artifacts.
- Full pipeline run with a mock backend.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple
from unittest.mock import MagicMock

import numpy as np
import pytest

from poltype.config.schema import PoltypeConfig
from poltype.errors import BackendError
from poltype.molecule.molecule import Molecule
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.factory import build_default_pipeline
from poltype.pipeline.runner import PipelineRunner
from poltype.pipeline.stage import StageResult, StageStatus
from poltype.pipeline.stages.esp import ESPFittingStage
from poltype.pipeline.stages.finalization import FinalizationStage
from poltype.pipeline.stages.geometry_opt import GeometryOptimizationStage
from poltype.pipeline.stages.multipole import MultipoleStage
from poltype.pipeline.stages.torsion import TorsionFittingStage
from poltype.qm.backend import (
    ESPGridResult,
    OptimizationResult,
    QMBackend,
    TorsionScanResult,
    WBOResult,
)


# ---------------------------------------------------------------------------
# Mock backend
# ---------------------------------------------------------------------------


class MockBackend(QMBackend):
    """In-memory mock backend that returns synthetic results."""

    def __init__(self) -> None:
        super().__init__(work_dir=None)
        self._opt_call_count = 0
        self._esp_call_count = 0
        self._torsion_call_count = 0
        self._wbo_call_count = 0

    @property
    def name(self) -> str:
        return "mock"

    def optimize_geometry(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> OptimizationResult:
        self._opt_call_count += 1
        return OptimizationResult(
            coordinates=molecule.coordinates,
            energy=-76.123456,
            converged=True,
        )

    def compute_esp_grid(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "aug-cc-pVTZ",
        **kwargs,
    ) -> ESPGridResult:
        self._esp_call_count += 1
        n_points = 50
        return ESPGridResult(
            grid_coords=np.random.randn(n_points, 3),
            esp_values=np.random.randn(n_points),
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
        self._torsion_call_count += 1
        scan_angles = angles if angles is not None else np.arange(-180, 180, 15.0)
        return TorsionScanResult(
            angles=scan_angles,
            energies=np.zeros_like(scan_angles),
        )

    def compute_wbo_matrix(
        self,
        molecule: Molecule,
        method: str = "HF",
        basis_set: str = "6-31G*",
        **kwargs,
    ) -> WBOResult:
        self._wbo_call_count += 1
        n = molecule.num_atoms
        return WBOResult(wbo_matrix=np.eye(n))


class FailingBackend(MockBackend):
    """Backend whose optimize_geometry raises a generic exception."""

    def optimize_geometry(self, molecule, method="MP2", basis_set="6-31G*", **kw):
        raise RuntimeError("QM calculation crashed")


class NonConvergingBackend(MockBackend):
    """Backend that returns a non-converged optimisation result."""

    def optimize_geometry(self, molecule, method="MP2", basis_set="6-31G*", **kw):
        return OptimizationResult(
            coordinates=molecule.coordinates,
            energy=-50.0,
            converged=False,
        )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def methane():
    return Molecule.from_smiles("C", name="methane", embed_3d=True)


@pytest.fixture
def ethanol():
    return Molecule.from_smiles("CCO", name="ethanol", embed_3d=True)


@pytest.fixture
def mock_backend():
    return MockBackend()


@pytest.fixture
def default_config():
    return PoltypeConfig()


@pytest.fixture
def context_with_backend(default_config, methane, mock_backend):
    return PipelineContext(
        config=default_config, molecule=methane, backend=mock_backend,
    )


@pytest.fixture
def ethanol_context(default_config, ethanol, mock_backend):
    return PipelineContext(
        config=default_config, molecule=ethanol, backend=mock_backend,
    )


# ===================================================================
# A. GeometryOptimizationStage tests
# ===================================================================


class TestGeometryOptimizationExecute:
    def test_successful_optimization(self, context_with_backend, mock_backend):
        stage = GeometryOptimizationStage()
        result = stage.execute(context_with_backend)
        assert result.status is StageStatus.COMPLETED
        assert "opt_result" in result.artifacts
        assert "molecule" in result.artifacts
        assert mock_backend._opt_call_count == 1

    def test_updates_molecule_in_context(self, context_with_backend):
        old_mol = context_with_backend.molecule
        stage = GeometryOptimizationStage()
        stage.execute(context_with_backend)
        # Molecule should have been replaced
        assert context_with_backend.molecule is not old_mol

    def test_energy_in_message(self, context_with_backend):
        stage = GeometryOptimizationStage()
        result = stage.execute(context_with_backend)
        assert "-76.123456" in result.message

    def test_non_converged_returns_failed(self, default_config, methane):
        backend = NonConvergingBackend()
        ctx = PipelineContext(
            config=default_config, molecule=methane, backend=backend,
        )
        stage = GeometryOptimizationStage()
        result = stage.execute(ctx)
        assert result.status is StageStatus.FAILED
        assert "did not converge" in result.message
        assert "opt_result" in result.artifacts

    def test_backend_error_raised(self, default_config, methane):
        backend = FailingBackend()
        ctx = PipelineContext(
            config=default_config, molecule=methane, backend=backend,
        )
        stage = GeometryOptimizationStage()
        with pytest.raises(BackendError, match="Geometry optimisation failed"):
            stage.execute(ctx)

    def test_not_implemented_propagates(self, default_config, methane):
        from poltype.qm.psi4_backend import Psi4Backend

        backend = Psi4Backend()
        ctx = PipelineContext(
            config=default_config, molecule=methane, backend=backend,
        )
        stage = GeometryOptimizationStage()
        with pytest.raises(NotImplementedError):
            stage.execute(ctx)

    def test_fails_when_no_molecule(self, default_config, mock_backend):
        ctx = PipelineContext(
            config=default_config, backend=mock_backend,
        )
        stage = GeometryOptimizationStage()
        result = stage.execute(ctx)
        assert result.status is StageStatus.FAILED

    def test_fails_when_no_backend(self, default_config, methane):
        ctx = PipelineContext(config=default_config, molecule=methane)
        stage = GeometryOptimizationStage()
        result = stage.execute(ctx)
        assert result.status is StageStatus.FAILED


# ===================================================================
# B. ESPFittingStage tests
# ===================================================================


class TestESPFittingExecute:
    def test_successful_esp(self, context_with_backend, mock_backend):
        stage = ESPFittingStage()
        result = stage.execute(context_with_backend)
        assert result.status is StageStatus.COMPLETED
        assert "esp_result" in result.artifacts
        assert mock_backend._esp_call_count == 1

    def test_esp_result_has_grid_points(self, context_with_backend):
        stage = ESPFittingStage()
        result = stage.execute(context_with_backend)
        esp = result.artifacts["esp_result"]
        assert len(esp.esp_values) == 50

    def test_message_contains_point_count(self, context_with_backend):
        stage = ESPFittingStage()
        result = stage.execute(context_with_backend)
        assert "50 points" in result.message

    def test_backend_error_raised_on_failure(self, default_config, methane):
        backend = MockBackend()
        backend.compute_esp_grid = MagicMock(
            side_effect=RuntimeError("ESP crashed")
        )
        ctx = PipelineContext(
            config=default_config, molecule=methane, backend=backend,
        )
        stage = ESPFittingStage()
        with pytest.raises(BackendError, match="ESP grid computation failed"):
            stage.execute(ctx)

    def test_not_implemented_propagates(self, default_config, methane):
        from poltype.qm.psi4_backend import Psi4Backend

        backend = Psi4Backend()
        ctx = PipelineContext(
            config=default_config, molecule=methane, backend=backend,
        )
        stage = ESPFittingStage()
        with pytest.raises(NotImplementedError):
            stage.execute(ctx)


# ===================================================================
# C. MultipoleStage tests
# ===================================================================


class TestMultipoleExecute:
    def test_completes_with_esp_artifact(self, context_with_backend):
        # First run ESP to produce the artifact
        esp_stage = ESPFittingStage()
        esp_result = esp_stage.execute(context_with_backend)
        context_with_backend.artifacts.update(esp_result.artifacts)

        stage = MultipoleStage()
        result = stage.execute(context_with_backend)
        assert result.status is StageStatus.COMPLETED
        assert "multipole_frames" in result.artifacts

    def test_multipole_frames_content(self, context_with_backend):
        esp_stage = ESPFittingStage()
        esp_result = esp_stage.execute(context_with_backend)
        context_with_backend.artifacts.update(esp_result.artifacts)

        stage = MultipoleStage()
        result = stage.execute(context_with_backend)
        frames = result.artifacts["multipole_frames"]
        assert frames["source"] == "esp"
        assert frames["num_atoms"] == context_with_backend.molecule.num_atoms

    def test_completes_without_esp_data(self, context_with_backend):
        stage = MultipoleStage()
        result = stage.execute(context_with_backend)
        assert result.status is StageStatus.COMPLETED
        assert "No ESP data" in result.message
        assert "multipole_frames" not in result.artifacts

    def test_fails_when_no_molecule(self, default_config, mock_backend):
        ctx = PipelineContext(config=default_config, backend=mock_backend)
        stage = MultipoleStage()
        result = stage.execute(ctx)
        assert result.status is StageStatus.FAILED


# ===================================================================
# D. TorsionFittingStage tests
# ===================================================================


class TestTorsionFittingExecute:
    def test_no_rotatable_bonds(self, context_with_backend):
        """Methane has no rotatable bonds."""
        stage = TorsionFittingStage()
        result = stage.execute(context_with_backend)
        assert result.status is StageStatus.COMPLETED
        assert "No rotatable bonds" in result.message
        assert result.artifacts["torsion_scans"] == {}

    def test_scans_rotatable_bonds(self, ethanol_context, mock_backend):
        """Ethanol has rotatable bonds; each should be scanned."""
        stage = TorsionFittingStage()
        result = stage.execute(ethanol_context)
        assert result.status is StageStatus.COMPLETED
        scans = result.artifacts["torsion_scans"]
        n_rotatable = len(ethanol_context.molecule.rotatable_bonds)
        assert len(scans) == n_rotatable
        assert mock_backend._torsion_call_count == n_rotatable

    def test_message_reports_count(self, ethanol_context):
        stage = TorsionFittingStage()
        result = stage.execute(ethanol_context)
        n = len(ethanol_context.molecule.rotatable_bonds)
        assert f"{n}" in result.message

    def test_backend_error_raised_on_failure(self, default_config, ethanol):
        backend = MockBackend()
        backend.torsion_scan = MagicMock(
            side_effect=RuntimeError("Torsion scan crashed")
        )
        ctx = PipelineContext(
            config=default_config, molecule=ethanol, backend=backend,
        )
        stage = TorsionFittingStage()
        with pytest.raises(BackendError, match="Torsion scan failed"):
            stage.execute(ctx)

    def test_not_implemented_propagates(self, default_config, ethanol):
        from poltype.qm.psi4_backend import Psi4Backend

        backend = Psi4Backend()
        ctx = PipelineContext(
            config=default_config, molecule=ethanol, backend=backend,
        )
        stage = TorsionFittingStage()
        with pytest.raises(NotImplementedError):
            stage.execute(ctx)


# ===================================================================
# E. FinalizationStage tests
# ===================================================================


class TestFinalizationExecute:
    def test_empty_context(self, default_config):
        ctx = PipelineContext(config=default_config)
        stage = FinalizationStage()
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        summary = result.artifacts["final_summary"]
        assert summary["has_molecule"] is False
        assert summary["collected_artifacts"] == []

    def test_collects_opt_result(self, context_with_backend):
        opt_result = OptimizationResult(
            coordinates=np.zeros((5, 3)), energy=-100.0,
        )
        context_with_backend.set_artifact("opt_result", opt_result)
        stage = FinalizationStage()
        result = stage.execute(context_with_backend)
        summary = result.artifacts["final_summary"]
        assert "opt_result" in summary["collected_artifacts"]
        assert summary["opt_energy"] == -100.0

    def test_collects_esp_result(self, context_with_backend):
        esp_result = ESPGridResult(
            grid_coords=np.zeros((25, 3)),
            esp_values=np.zeros(25),
        )
        context_with_backend.set_artifact("esp_result", esp_result)
        stage = FinalizationStage()
        result = stage.execute(context_with_backend)
        summary = result.artifacts["final_summary"]
        assert "esp_result" in summary["collected_artifacts"]
        assert summary["esp_grid_points"] == 25

    def test_collects_torsion_scans(self, context_with_backend):
        context_with_backend.set_artifact("torsion_scans", {(0, 1): "scan"})
        stage = FinalizationStage()
        result = stage.execute(context_with_backend)
        summary = result.artifacts["final_summary"]
        assert "torsion_scans" in summary["collected_artifacts"]
        assert summary["torsion_bonds_scanned"] == 1

    def test_collects_multipole_frames(self, context_with_backend):
        frames = {"source": "esp", "num_atoms": 5}
        context_with_backend.set_artifact("multipole_frames", frames)
        stage = FinalizationStage()
        result = stage.execute(context_with_backend)
        summary = result.artifacts["final_summary"]
        assert "multipole_frames" in summary["collected_artifacts"]

    def test_molecule_info_in_summary(self, context_with_backend):
        stage = FinalizationStage()
        result = stage.execute(context_with_backend)
        summary = result.artifacts["final_summary"]
        assert summary["molecule_name"] == "methane"
        assert summary["num_atoms"] > 0
        assert summary["formula"] != ""


# ===================================================================
# F. Full pipeline with mock backend
# ===================================================================


class TestFullPipelineWithMockBackend:
    def test_full_run_completes(self, default_config, methane, mock_backend):
        runner = build_default_pipeline()
        ctx = PipelineContext(
            config=default_config,
            molecule=methane,
            backend=mock_backend,
        )
        result = runner.run(ctx)
        # All 6 stages should complete
        assert len(result.stage_results) == 6
        for sr in result.stage_results.values():
            assert sr.status is StageStatus.COMPLETED

    def test_final_summary_present(self, default_config, methane, mock_backend):
        runner = build_default_pipeline()
        ctx = PipelineContext(
            config=default_config,
            molecule=methane,
            backend=mock_backend,
        )
        result = runner.run(ctx)
        summary = result.get_artifact("final_summary")
        assert summary is not None
        assert summary["has_molecule"] is True
        assert "opt_result" in summary["collected_artifacts"]
        assert "esp_result" in summary["collected_artifacts"]

    def test_database_match_only_skips_qm_stages(self, default_config, methane):
        default_config.database_match_only = True
        runner = build_default_pipeline()
        ctx = PipelineContext(config=default_config, molecule=methane)
        result = runner.run(ctx)
        # input_preparation and finalization should complete
        assert (
            result.stage_results["input_preparation"].status
            is StageStatus.COMPLETED
        )
        # QM stages should be skipped
        assert (
            result.stage_results["geometry_optimization"].status
            is StageStatus.SKIPPED
        )
        assert (
            result.stage_results["esp_fitting"].status is StageStatus.SKIPPED
        )
        assert (
            result.stage_results["multipole"].status is StageStatus.SKIPPED
        )
        assert (
            result.stage_results["torsion_fitting"].status
            is StageStatus.SKIPPED
        )
        assert (
            result.stage_results["finalization"].status
            is StageStatus.COMPLETED
        )

    def test_opt_energy_propagates_to_finalization(
        self, default_config, methane, mock_backend,
    ):
        runner = build_default_pipeline()
        ctx = PipelineContext(
            config=default_config,
            molecule=methane,
            backend=mock_backend,
        )
        result = runner.run(ctx)
        summary = result.get_artifact("final_summary")
        assert summary["opt_energy"] == pytest.approx(-76.123456)

    def test_ethanol_torsion_scans_collected(
        self, default_config, ethanol, mock_backend,
    ):
        runner = build_default_pipeline()
        ctx = PipelineContext(
            config=default_config,
            molecule=ethanol,
            backend=mock_backend,
        )
        result = runner.run(ctx)
        summary = result.get_artifact("final_summary")
        assert summary["torsion_bonds_scanned"] >= 1

    def test_non_converged_opt_halts_pipeline(self, default_config, methane):
        backend = NonConvergingBackend()
        runner = build_default_pipeline()
        ctx = PipelineContext(
            config=default_config,
            molecule=methane,
            backend=backend,
        )
        result = runner.run(ctx)
        assert (
            result.stage_results["geometry_optimization"].status
            is StageStatus.FAILED
        )
        # Pipeline should stop — no esp_fitting result
        assert "esp_fitting" not in result.stage_results
