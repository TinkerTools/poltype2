"""
tests/unit/test_phase4.py – unit tests for Phase 4: stage implementations.

Tests cover:
- GeometryOptimizationStage: delegates to backend, updates molecule,
  handles convergence failure and backend errors.
- MultipoleStage: runs DMA to produce initial multipole frames.
- ESPFittingStage: delegates to backend, stores ESP artifact, optimises
  multipoles from the upstream MultipoleStage.
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
    DMAResult,
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

    def compute_dma(
        self,
        molecule: Molecule,
        method: str = "MP2",
        basis_set: str = "6-311G**",
        **kwargs,
    ) -> DMAResult:
        from pathlib import Path
        fchk = Path(molecule.work_dir) / f"{molecule.name or 'mol'}_dma.fchk"
        fchk.write_text("mock fchk")
        return DMAResult(fchk_path=fchk, energy=-76.0)

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

    def test_generates_fitted_parameter_artifacts(self, context_with_backend):
        stage = ESPFittingStage()
        result = stage.execute(context_with_backend)
        assert "fitted_multipoles_by_atom_index" in result.artifacts
        assert "polarization_params_by_atom_index" in result.artifacts
        assert len(result.artifacts["fitted_multipoles_by_atom_index"]) == (
            context_with_backend.molecule.num_atoms
        )

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
    def test_completes_with_backend(self, context_with_backend):
        stage = MultipoleStage()
        result = stage.execute(context_with_backend)
        assert result.status is StageStatus.COMPLETED
        assert "multipole_frames" in result.artifacts

    def test_multipole_frames_content(self, context_with_backend):
        stage = MultipoleStage()
        result = stage.execute(context_with_backend)
        frames = result.artifacts["multipole_frames"]
        assert frames["source"] == "dma"
        assert frames["num_atoms"] == context_with_backend.molecule.num_atoms
        assert frames["dma_method"] == "MP2"
        assert frames["dma_basis_set"] == "6-311G**"

    def test_fails_when_no_molecule(self, default_config, mock_backend):
        ctx = PipelineContext(config=default_config, backend=mock_backend)
        stage = MultipoleStage()
        result = stage.execute(ctx)
        assert result.status is StageStatus.FAILED

    def test_fails_when_no_backend(self, default_config, methane):
        ctx = PipelineContext(config=default_config, molecule=methane)
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

    def test_skips_database_matched_torsions(self, ethanol_context, mock_backend):
        matched_torsions = {
            (0, bond[0], bond[1], 0)
            for bond in ethanol_context.molecule.rotatable_bonds
        }
        ethanol_context.set_artifact("db_match_result", MagicMock(
            matched_torsions=matched_torsions
        ))
        stage = TorsionFittingStage()
        result = stage.execute(ethanol_context)
        assert result.status is StageStatus.COMPLETED
        assert result.artifacts["torsion_scans"] == {}
        assert result.artifacts["torsion_params"] == []
        assert mock_backend._torsion_call_count == 0

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
        frames = {"source": "dma", "num_atoms": 5}
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
        # All 10 stages should complete
        assert len(result.stage_results) == 10
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
        assert summary["torsion_bonds_scanned"] >= 0

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
        # Pipeline should stop — no multipole result
        assert "multipole" not in result.stage_results


# ===================================================================
# G. New workflow integration tests (Steps 0-6)
# ===================================================================


class TestOptXyzWriting:
    """Tests for writing {fname}_opt.xyz after geometry optimisation."""

    def test_opt_xyz_written(self, context_with_backend, tmp_path):
        """Geometry optimisation should write {fname}_opt.xyz."""
        context_with_backend._work_dir = tmp_path
        stage = GeometryOptimizationStage()
        result = stage.execute(context_with_backend)
        assert result.status is StageStatus.COMPLETED
        assert "opt_xyz_path" in result.artifacts
        opt_xyz = result.artifacts["opt_xyz_path"]
        assert opt_xyz.exists()
        assert opt_xyz.name == "methane_opt.xyz"

    def test_opt_xyz_content(self, context_with_backend, tmp_path):
        """The opt.xyz file should contain valid atom lines."""
        context_with_backend._work_dir = tmp_path
        stage = GeometryOptimizationStage()
        result = stage.execute(context_with_backend)
        text = result.artifacts["opt_xyz_path"].read_text()
        lines = text.strip().splitlines()
        # First line is the header with atom count
        assert "5" in lines[0]  # methane has 5 atoms
        assert "methane" in lines[0]
        # Data lines should have atom index, symbol, coords
        assert len(lines) == 6  # 1 header + 5 atoms

    def test_opt_xyz_coordinates_match(self, context_with_backend, tmp_path):
        """Coordinates in opt.xyz should match the optimised molecule."""
        context_with_backend._work_dir = tmp_path
        stage = GeometryOptimizationStage()
        result = stage.execute(context_with_backend)
        mol = result.artifacts["molecule"]
        text = result.artifacts["opt_xyz_path"].read_text()
        lines = text.strip().splitlines()
        # Parse first atom coordinates from line 2
        parts = lines[1].split()
        # parts: index, symbol, x, y, z, type, neighbor(s)
        x_file = float(parts[2])
        x_mol = mol.coordinates[0][0]
        assert abs(x_file - x_mol) < 1e-4


class TestDMAResultInMultipoleStage:
    """Tests for DMA QM calculation in MultipoleStage."""

    def test_dma_result_artifact_present(self, context_with_backend):
        """MultipoleStage should produce a dma_result artifact."""
        stage = MultipoleStage()
        result = stage.execute(context_with_backend)
        assert result.status is StageStatus.COMPLETED
        assert "dma_result" in result.artifacts

    def test_dma_result_has_fchk_path(self, context_with_backend):
        """The DMA result should contain a path to the fchk file."""
        stage = MultipoleStage()
        result = stage.execute(context_with_backend)
        dma_result = result.artifacts["dma_result"]
        assert dma_result.fchk_path is not None
        assert dma_result.fchk_path.exists()

    def test_dma_result_energy(self, context_with_backend):
        """The DMA result should contain an energy value."""
        stage = MultipoleStage()
        result = stage.execute(context_with_backend)
        dma_result = result.artifacts["dma_result"]
        assert dma_result.energy == pytest.approx(-76.0)

    def test_multipole_frames_contain_fchk_path(self, context_with_backend):
        """The multipole_frames artifact should record the fchk path."""
        stage = MultipoleStage()
        result = stage.execute(context_with_backend)
        frames = result.artifacts["multipole_frames"]
        assert "fchk_path" in frames


class TestESPGridWriting:
    """Tests for ESP grid file writing in ESPFittingStage."""

    def test_esp_grid_file_written(self, context_with_backend, tmp_path):
        """ESPFittingStage should write an ESP grid file."""
        context_with_backend._work_dir = tmp_path
        stage = ESPFittingStage()
        result = stage.execute(context_with_backend)
        assert result.status is StageStatus.COMPLETED
        assert "esp_grid_path" in result.artifacts
        grid_path = result.artifacts["esp_grid_path"]
        assert grid_path.exists()
        assert grid_path.name == "methane_esp.grid"

    def test_esp_grid_content(self, context_with_backend, tmp_path):
        """The ESP grid file should contain the correct number of points."""
        context_with_backend._work_dir = tmp_path
        stage = ESPFittingStage()
        result = stage.execute(context_with_backend)
        grid_path = result.artifacts["esp_grid_path"]
        text = grid_path.read_text()
        lines = text.strip().splitlines()
        # First line is the point count
        n_points = int(lines[0].strip())
        assert n_points == 50
        # Subsequent lines have x, y, z, potential
        assert len(lines) == n_points + 1


class TestDMAResultDataclass:
    """Tests for the DMAResult dataclass itself."""

    def test_construction(self, tmp_path):
        """DMAResult can be constructed with required fields."""
        fchk = tmp_path / "test.fchk"
        fchk.write_text("dummy")
        result = DMAResult(fchk_path=fchk, energy=-100.0)
        assert result.fchk_path == fchk
        assert result.energy == -100.0

    def test_default_log_path_is_none(self, tmp_path):
        """DMAResult.log_path defaults to None."""
        fchk = tmp_path / "test.fchk"
        fchk.write_text("dummy")
        result = DMAResult(fchk_path=fchk)
        assert result.log_path is None

    def test_default_energy_is_zero(self, tmp_path):
        """DMAResult.energy defaults to 0.0."""
        fchk = tmp_path / "test.fchk"
        fchk.write_text("dummy")
        result = DMAResult(fchk_path=fchk)
        assert result.energy == 0.0


class TestExternalToolRunners:
    """Tests for external tool runner construction and input generation."""

    def test_gdma_runner_construction(self):
        from poltype.external.gdma import GDMARunner
        runner = GDMARunner(gdma_exe="/usr/local/bin/gdma")
        assert runner.gdma_exe == "/usr/local/bin/gdma"

    def test_gdma_input_generation(self, tmp_path):
        from poltype.external.gdma import GDMARunner
        fchk_path = tmp_path / "test.fchk"
        punch_path = tmp_path / "test.punch"
        inp = GDMARunner._build_input(fchk_path, punch_path, multipole_rank=2)
        assert "File" in inp
        assert fchk_path.name in inp
        assert str(fchk_path) not in inp
        assert "Punch" in inp
        assert punch_path.name in inp
        assert str(punch_path) not in inp
        assert "Limit 2" in inp

    def test_poledit_runner_construction(self):
        from poltype.external.poledit import PoleditRunner
        runner = PoleditRunner(poledit_exe="/usr/local/bin/poledit")
        assert runner.poledit_exe == "/usr/local/bin/poledit"

    def test_potential_runner_construction(self):
        from poltype.external.potential import PotentialRunner
        runner = PotentialRunner(potential_exe="/usr/local/bin/potential")
        assert runner.potential_exe == "/usr/local/bin/potential"

    def test_gdma_result_dataclass(self, tmp_path):
        from poltype.external.gdma import GDMAResult
        punch = tmp_path / "test.punch"
        punch.write_text("dummy")
        result = GDMAResult(punch_path=punch)
        assert result.punch_path == punch
        assert result.log_path is None

    def test_gdma_rejects_missing_fchk(self, tmp_path):
        from poltype.external.gdma import GDMARunner
        runner = GDMARunner()
        with pytest.raises(FileNotFoundError, match="does not exist"):
            runner.run(
                fchk_path=tmp_path / "missing.fchk",
                work_dir=tmp_path,
            )

    def test_gdma_rejects_non_fchk_extension(self, tmp_path):
        from poltype.external.gdma import GDMARunner
        molden = tmp_path / "dma.molden"
        molden.write_text("dummy")
        runner = GDMARunner()
        with pytest.raises(ValueError, match="requires a .fchk file"):
            runner.run(fchk_path=molden, work_dir=tmp_path)

    def test_poledit_result_dataclass(self):
        from poltype.external.poledit import PoleditResult
        result = PoleditResult()
        assert result.key_path is None
        assert result.log_path is None

    def test_potential_result_dataclass(self):
        from poltype.external.potential import PotentialResult
        result = PotentialResult()
        assert result.key_path is None
        assert result.log_path is None


class TestFullPipelineWorkflow:
    """Integration tests for the full poltype workflow steps 0-6."""

    def test_opt_xyz_available_in_pipeline(self, default_config, methane, mock_backend, tmp_path):
        """Full pipeline should produce opt_xyz_path artifact."""
        runner = build_default_pipeline()
        ctx = PipelineContext(
            config=default_config,
            molecule=methane,
            backend=mock_backend,
            work_dir=tmp_path,
        )
        result = runner.run(ctx)
        assert "opt_xyz_path" in result.artifacts
        assert result.artifacts["opt_xyz_path"].exists()

    def test_dma_result_available_in_pipeline(self, default_config, methane, mock_backend, tmp_path):
        """Full pipeline should produce dma_result artifact."""
        runner = build_default_pipeline()
        ctx = PipelineContext(
            config=default_config,
            molecule=methane,
            backend=mock_backend,
            work_dir=tmp_path,
        )
        result = runner.run(ctx)
        assert "dma_result" in result.artifacts

    def test_esp_grid_available_in_pipeline(self, default_config, methane, mock_backend, tmp_path):
        """Full pipeline should produce esp_grid_path artifact."""
        runner = build_default_pipeline()
        ctx = PipelineContext(
            config=default_config,
            molecule=methane,
            backend=mock_backend,
            work_dir=tmp_path,
        )
        result = runner.run(ctx)
        assert "esp_grid_path" in result.artifacts
        assert result.artifacts["esp_grid_path"].exists()
