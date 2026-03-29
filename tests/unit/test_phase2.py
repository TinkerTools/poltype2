"""
tests/unit/test_phase2.py – unit tests for Phase 2: Pipeline framework
and QM backend abstraction.

Tests cover:
- StageStatus enum values and member count.
- StageResult dataclass construction, defaults, and artifact isolation.
- Stage ABC enforcement and concrete subclass behaviour.
- PipelineContext construction, artifact helpers, molecule updates,
  work_dir handling, and repr.
- PipelineRunner ordering, failure halting, skip handling, artifact
  merging, elapsed-time tracking, and chaining API.
- Concrete stage stub names, execute() preconditions, and should_skip()
  logic for InputPreparation, GeometryOptimization, ESPFitting,
  Multipole, TorsionFitting, and Finalization stages.
- QMBackend ABC enforcement and partial-implementation guard.
- Result dataclasses (OptimizationResult, ESPGridResult,
  TorsionScanResult, WBOResult).
- Concrete backend stubs (Psi4, Gaussian, PySCF) and their interface.
- select_backend() factory: config-driven backend selection, argument
  propagation, and repr.
- Pipeline ↔ Backend integration: context carries backend, stages
  interact with backend, and backend selection from config.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from poltype.config.schema import PoltypeConfig, QMConfig, ResourceConfig
from poltype.molecule.molecule import Molecule
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.runner import PipelineRunner
from poltype.pipeline.stage import Stage, StageResult, StageStatus
from poltype.pipeline.stages.esp import ESPFittingStage
from poltype.pipeline.stages.finalization import FinalizationStage
from poltype.pipeline.stages.geometry_opt import GeometryOptimizationStage
from poltype.pipeline.stages.input_prep import InputPreparationStage
from poltype.pipeline.stages.multipole import MultipoleStage
from poltype.pipeline.stages.torsion import TorsionFittingStage
from poltype.qm.backend import (
    ESPGridResult,
    OptimizationResult,
    QMBackend,
    TorsionScanResult,
    WBOResult,
)
from poltype.qm.factory import select_backend
from poltype.qm.gaussian_backend import GaussianBackend
from poltype.qm.psi4_backend import Psi4Backend
from poltype.qm.pyscf_backend import PySCFBackend


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def default_config():
    return PoltypeConfig()


@pytest.fixture
def methane():
    return Molecule.from_smiles("C", name="methane", embed_3d=True)


@pytest.fixture
def context(default_config):
    return PipelineContext(config=default_config)


@pytest.fixture
def context_with_molecule(default_config, methane):
    return PipelineContext(config=default_config, molecule=methane)


# ---------------------------------------------------------------------------
# Helpers — simple concrete stages for testing the runner
# ---------------------------------------------------------------------------


class _CountingStage(Stage):
    """Stage that appends its name to an 'order' list in artifacts."""

    def execute(self, context: PipelineContext) -> StageResult:
        order = context.get_artifact("order", [])
        order.append(self.name)
        return StageResult(
            status=StageStatus.COMPLETED,
            artifacts={"order": order},
        )


class _FailingStage(Stage):
    """Stage that always fails."""

    def execute(self, context: PipelineContext) -> StageResult:
        return StageResult(
            status=StageStatus.FAILED,
            message="intentional failure",
        )


class _AlwaysSkipStage(Stage):
    """Stage whose should_skip always returns True."""

    def should_skip(self, context: PipelineContext) -> bool:
        return True

    def execute(self, context: PipelineContext) -> StageResult:
        return StageResult(
            status=StageStatus.COMPLETED,
            artifacts={"should_not_appear": True},
        )


class _ArtifactStage(Stage):
    """Stage that produces a specific artifact."""

    def __init__(self, name: str, key: str, value: object) -> None:
        super().__init__(name)
        self._key = key
        self._value = value

    def execute(self, context: PipelineContext) -> StageResult:
        return StageResult(
            status=StageStatus.COMPLETED,
            artifacts={self._key: self._value},
        )


# ===================================================================
# A. StageStatus enum
# ===================================================================


class TestStageStatus:
    def test_all_enum_values(self):
        assert StageStatus.PENDING.value == "pending"
        assert StageStatus.RUNNING.value == "running"
        assert StageStatus.COMPLETED.value == "completed"
        assert StageStatus.SKIPPED.value == "skipped"
        assert StageStatus.FAILED.value == "failed"

    def test_enum_member_count(self):
        assert len(StageStatus) == 5


# ===================================================================
# B. StageResult dataclass
# ===================================================================


class TestStageResult:
    def test_status_only(self):
        r = StageResult(status=StageStatus.COMPLETED)
        assert r.status is StageStatus.COMPLETED
        assert r.message == ""
        assert r.elapsed_seconds == 0.0
        assert r.artifacts == {}

    def test_all_fields(self):
        r = StageResult(
            status=StageStatus.FAILED,
            message="boom",
            elapsed_seconds=1.23,
            artifacts={"k": "v"},
        )
        assert r.status is StageStatus.FAILED
        assert r.message == "boom"
        assert r.elapsed_seconds == 1.23
        assert r.artifacts == {"k": "v"}

    def test_artifacts_default_independent(self):
        r1 = StageResult(status=StageStatus.COMPLETED)
        r2 = StageResult(status=StageStatus.COMPLETED)
        r1.artifacts["a"] = 1
        assert "a" not in r2.artifacts


# ===================================================================
# C. Stage ABC
# ===================================================================


class TestStageABC:
    def test_cannot_instantiate_directly(self):
        with pytest.raises(TypeError):
            Stage("some-stage")  # type: ignore[abstract]

    def test_concrete_subclass(self):
        stage = _CountingStage("counter")
        assert stage.name == "counter"

    def test_should_skip_default_false(self, context):
        stage = _CountingStage("counter")
        assert stage.should_skip(context) is False

    def test_name_property(self):
        stage = _CountingStage("my-stage")
        assert stage.name == "my-stage"

    def test_repr(self):
        stage = _CountingStage("repr-test")
        assert repr(stage) == "Stage(name='repr-test')"


# ===================================================================
# D. PipelineContext
# ===================================================================


class TestPipelineContext:
    def test_config_only(self, default_config):
        ctx = PipelineContext(config=default_config)
        assert ctx.config is default_config
        assert ctx.molecule is None
        assert ctx.backend is None
        assert ctx.work_dir == Path.cwd()

    def test_molecule_starts_none(self, default_config):
        ctx = PipelineContext(config=default_config)
        assert ctx.molecule is None

    def test_work_dir_defaults_to_cwd(self, default_config):
        ctx = PipelineContext(config=default_config)
        assert ctx.work_dir == Path.cwd()

    def test_custom_work_dir(self, default_config, tmp_path):
        ctx = PipelineContext(config=default_config, work_dir=tmp_path)
        assert ctx.work_dir == tmp_path.resolve()

    def test_set_and_get_artifact(self, context):
        context.set_artifact("opt_result", {"energy": -100.0})
        assert context.get_artifact("opt_result") == {"energy": -100.0}

    def test_get_artifact_missing_returns_default(self, context):
        assert context.get_artifact("missing") is None
        assert context.get_artifact("missing", 42) == 42

    def test_update_molecule(self, context, methane):
        assert context.molecule is None
        context.update_molecule(methane)
        assert context.molecule is methane

    def test_stage_results_starts_empty(self, context):
        assert context.stage_results == {}

    def test_molecule_setter(self, context, methane):
        context.molecule = methane
        assert context.molecule is methane

    def test_repr(self, context):
        r = repr(context)
        assert "PipelineContext" in r

    def test_repr_with_molecule(self, context_with_molecule):
        assert "methane" in repr(context_with_molecule)

    def test_backend_initially_none(self, context):
        assert context.backend is None

    def test_backend_setter(self, context, default_config):
        backend = select_backend(default_config)
        context.backend = backend
        assert context.backend is backend


# ===================================================================
# E. PipelineRunner
# ===================================================================


class TestPipelineRunner:
    def test_empty_runner(self):
        runner = PipelineRunner()
        assert repr(runner) == "PipelineRunner(stages=0)"

    def test_add_stage_returns_self(self):
        runner = PipelineRunner()
        ret = runner.add_stage(_CountingStage("a"))
        assert ret is runner

    def test_chaining(self):
        runner = (
            PipelineRunner()
            .add_stage(_CountingStage("a"))
            .add_stage(_CountingStage("b"))
        )
        assert repr(runner) == "PipelineRunner(stages=2)"

    def test_run_empty(self, context):
        runner = PipelineRunner()
        result = runner.run(context)
        assert result is context
        assert result.stage_results == {}

    def test_run_order(self, context):
        runner = PipelineRunner(
            stages=[
                _CountingStage("first"),
                _CountingStage("second"),
                _CountingStage("third"),
            ]
        )
        result = runner.run(context)
        assert result.get_artifact("order") == ["first", "second", "third"]

    def test_run_records_results(self, context):
        runner = PipelineRunner(
            stages=[_CountingStage("a"), _CountingStage("b")]
        )
        result = runner.run(context)
        assert "a" in result.stage_results
        assert "b" in result.stage_results
        assert result.stage_results["a"].status is StageStatus.COMPLETED

    def test_run_stops_on_failure(self, context):
        runner = PipelineRunner(
            stages=[
                _CountingStage("before"),
                _FailingStage("boom"),
                _CountingStage("after"),
            ]
        )
        result = runner.run(context)
        assert "before" in result.stage_results
        assert "boom" in result.stage_results
        assert "after" not in result.stage_results
        assert result.stage_results["boom"].status is StageStatus.FAILED

    def test_run_skips_stages(self, context):
        runner = PipelineRunner(
            stages=[
                _CountingStage("before"),
                _AlwaysSkipStage("skipped"),
                _CountingStage("after"),
            ]
        )
        result = runner.run(context)
        assert result.stage_results["skipped"].status is StageStatus.SKIPPED
        assert "should_not_appear" not in result.artifacts
        assert result.get_artifact("order") == ["before", "after"]

    def test_run_merges_artifacts(self, context):
        runner = PipelineRunner(
            stages=[
                _ArtifactStage("s1", "key1", "val1"),
                _ArtifactStage("s2", "key2", "val2"),
            ]
        )
        result = runner.run(context)
        assert result.get_artifact("key1") == "val1"
        assert result.get_artifact("key2") == "val2"

    def test_elapsed_seconds_populated(self, context):
        runner = PipelineRunner(stages=[_CountingStage("timed")])
        result = runner.run(context)
        assert result.stage_results["timed"].elapsed_seconds >= 0.0

    def test_instantiate_with_stages_list(self):
        stages = [_CountingStage("a"), _CountingStage("b")]
        runner = PipelineRunner(stages=stages)
        assert repr(runner) == "PipelineRunner(stages=2)"


# ===================================================================
# F. Concrete stage stubs
# ===================================================================


class TestInputPreparationStage:
    def test_name(self):
        assert InputPreparationStage().name == "input_preparation"

    def test_execute_with_molecule(self, context_with_molecule):
        result = InputPreparationStage().execute(context_with_molecule)
        assert result.status is StageStatus.COMPLETED
        assert "already present" in result.message

    def test_execute_no_structure_fails(self, context):
        result = InputPreparationStage().execute(context)
        assert result.status is StageStatus.FAILED
        assert "No structure" in result.message

    def test_should_skip_check_input_only(self, default_config):
        default_config.check_input_only = True
        ctx = PipelineContext(config=default_config)
        assert InputPreparationStage().should_skip(ctx) is True

    def test_should_not_skip_default(self, context):
        assert InputPreparationStage().should_skip(context) is False


class TestGeometryOptimizationStage:
    def test_name(self):
        assert GeometryOptimizationStage().name == "geometry_optimization"

    def test_fails_no_molecule(self, context):
        result = GeometryOptimizationStage().execute(context)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_fails_no_backend(self, context_with_molecule):
        result = GeometryOptimizationStage().execute(context_with_molecule)
        assert result.status is StageStatus.FAILED
        assert "No QM backend" in result.message

    def test_skip_database_match_only(self, default_config):
        default_config.database_match_only = True
        ctx = PipelineContext(config=default_config)
        assert GeometryOptimizationStage().should_skip(ctx) is True

    def test_should_not_skip_default(self, context):
        assert GeometryOptimizationStage().should_skip(context) is False


class TestESPFittingStage:
    def test_name(self):
        assert ESPFittingStage().name == "esp_fitting"

    def test_fails_no_molecule(self, context):
        result = ESPFittingStage().execute(context)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_fails_no_backend(self, context_with_molecule):
        result = ESPFittingStage().execute(context_with_molecule)
        assert result.status is StageStatus.FAILED
        assert "No QM backend" in result.message

    def test_skip_database_match_only(self, default_config):
        default_config.database_match_only = True
        ctx = PipelineContext(config=default_config)
        assert ESPFittingStage().should_skip(ctx) is True

    def test_should_not_skip_default(self, context):
        assert ESPFittingStage().should_skip(context) is False


class TestMultipoleStage:
    def test_name(self):
        assert MultipoleStage().name == "multipole"

    def test_fails_no_molecule(self, context):
        result = MultipoleStage().execute(context)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_completes_with_molecule_and_backend(self, context_with_molecule):
        from unittest.mock import MagicMock
        context_with_molecule._backend = MagicMock()
        context_with_molecule._backend.name = "mock"
        result = MultipoleStage().execute(context_with_molecule)
        assert result.status is StageStatus.COMPLETED
        assert "multipole_frames" in result.artifacts

    def test_fails_without_backend(self, context_with_molecule):
        result = MultipoleStage().execute(context_with_molecule)
        assert result.status is StageStatus.FAILED
        assert "No QM backend" in result.message

    def test_skip_database_match_only(self, default_config):
        default_config.database_match_only = True
        ctx = PipelineContext(config=default_config)
        assert MultipoleStage().should_skip(ctx) is True

    def test_should_not_skip_default(self, context):
        assert MultipoleStage().should_skip(context) is False


class TestTorsionFittingStage:
    def test_name(self):
        assert TorsionFittingStage().name == "torsion_fitting"

    def test_fails_no_molecule(self, context):
        result = TorsionFittingStage().execute(context)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_fails_no_backend(self, context_with_molecule):
        result = TorsionFittingStage().execute(context_with_molecule)
        assert result.status is StageStatus.FAILED
        assert "No QM backend" in result.message

    def test_skip_torsion_fit_false(self, default_config):
        default_config.do_torsion_fit = False
        ctx = PipelineContext(config=default_config)
        assert TorsionFittingStage().should_skip(ctx) is True

    def test_skip_database_match_only(self, default_config):
        default_config.database_match_only = True
        ctx = PipelineContext(config=default_config)
        assert TorsionFittingStage().should_skip(ctx) is True

    def test_should_not_skip_default(self, context):
        assert TorsionFittingStage().should_skip(context) is False


class TestFinalizationStage:
    def test_name(self):
        assert FinalizationStage().name == "finalization"

    def test_execute_completed(self, context):
        result = FinalizationStage().execute(context)
        assert result.status is StageStatus.COMPLETED
        assert "final_summary" in result.artifacts

    def test_execute_collects_molecule_info(self, context_with_molecule):
        result = FinalizationStage().execute(context_with_molecule)
        summary = result.artifacts["final_summary"]
        assert summary["has_molecule"] is True
        assert summary["molecule_name"] == "methane"


# ===================================================================
# G. QMBackend ABC
# ===================================================================


class TestQMBackendABC:
    def test_cannot_instantiate(self):
        with pytest.raises(TypeError):
            QMBackend()  # type: ignore[abstract]

    def test_partial_implementation_fails(self):
        class PartialBackend(QMBackend):
            def optimize_geometry(self, *a, **kw):
                pass

        with pytest.raises(TypeError):
            PartialBackend()


# ===================================================================
# H. Result dataclasses
# ===================================================================


class TestResultDataclasses:
    def test_optimization_result(self):
        coords = np.zeros((5, 3))
        r = OptimizationResult(coordinates=coords, energy=-100.0)
        assert r.converged is True
        assert r.log_path is None

    def test_esp_grid_result(self):
        grid = np.zeros((100, 3))
        vals = np.zeros(100)
        r = ESPGridResult(grid_coords=grid, esp_values=vals)
        assert r.log_path is None

    def test_torsion_scan_result(self):
        angles = np.arange(-180, 180, 15, dtype=float)
        energies = np.zeros_like(angles)
        r = TorsionScanResult(angles=angles, energies=energies)
        assert len(r.log_paths) == 0

    def test_wbo_result(self):
        mat = np.eye(5)
        r = WBOResult(wbo_matrix=mat)
        assert r.log_path is None


# ===================================================================
# I. Concrete backend stubs
# ===================================================================


class TestPsi4Backend:
    def test_instantiation(self):
        assert isinstance(Psi4Backend(), QMBackend)

    def test_name(self):
        assert Psi4Backend().name == "psi4"

    def test_optimize_raises(self):
        mol = Molecule.from_smiles("C", name="m", embed_3d=True)
        with pytest.raises(NotImplementedError):
            Psi4Backend().optimize_geometry(mol, "MP2", "6-31G*")

    def test_esp_raises(self):
        mol = Molecule.from_smiles("C", name="m", embed_3d=True)
        with pytest.raises(NotImplementedError):
            Psi4Backend().compute_esp_grid(mol, "MP2", "aug-cc-pVTZ")

    def test_torsion_raises(self):
        mol = Molecule.from_smiles("C", name="m", embed_3d=True)
        with pytest.raises(NotImplementedError):
            Psi4Backend().torsion_scan(mol, (0, 1), "wB97X-D", "6-31G*")

    def test_wbo_raises(self):
        mol = Molecule.from_smiles("C", name="m", embed_3d=True)
        with pytest.raises(NotImplementedError):
            Psi4Backend().compute_wbo_matrix(mol, "HF", "6-31G*")


class TestGaussianBackend:
    def test_instantiation(self):
        assert isinstance(GaussianBackend(), QMBackend)

    def test_name(self):
        assert GaussianBackend().name == "gaussian"

    def test_optimize_raises(self):
        mol = Molecule.from_smiles("C", name="m", embed_3d=True)
        with pytest.raises(NotImplementedError):
            GaussianBackend().optimize_geometry(mol, "MP2", "6-31G*")


class TestPySCFBackend:
    def test_instantiation(self):
        assert isinstance(PySCFBackend(), QMBackend)

    def test_name(self):
        assert PySCFBackend().name == "pyscf"

    def test_optimize_raises(self):
        mol = Molecule.from_smiles("C", name="m", embed_3d=True)
        with pytest.raises(NotImplementedError):
            PySCFBackend().optimize_geometry(mol, "MP2", "6-31G*")


# ===================================================================
# J. select_backend factory
# ===================================================================


class TestSelectBackend:
    def test_default_is_pyscf(self, default_config):
        backend = select_backend(default_config)
        assert isinstance(backend, PySCFBackend)

    def test_dont_use_pyscf_selects_psi4(self):
        cfg = PoltypeConfig(qm=QMConfig(dont_use_pyscf=True))
        assert isinstance(select_backend(cfg), Psi4Backend)

    def test_use_gaus_selects_gaussian(self):
        cfg = PoltypeConfig(qm=QMConfig(use_gaus=True))
        assert isinstance(select_backend(cfg), GaussianBackend)

    def test_use_gaus_opt_only_selects_gaussian(self):
        cfg = PoltypeConfig(qm=QMConfig(use_gaus_opt_only=True))
        assert isinstance(select_backend(cfg), GaussianBackend)

    def test_psi4_args_propagated(self):
        cfg = PoltypeConfig(
            qm=QMConfig(dont_use_pyscf=True, psi4_args="--test")
        )
        backend = select_backend(cfg)
        assert isinstance(backend, Psi4Backend)
        assert backend.psi4_args == "--test"

    def test_num_proc_propagated(self):
        cfg = PoltypeConfig(
            qm=QMConfig(dont_use_pyscf=True),
            resources=ResourceConfig(num_proc=8),
        )
        backend = select_backend(cfg)
        assert backend.num_proc == 8

    def test_repr_contains_name(self, default_config):
        backend = select_backend(default_config)
        assert backend.name in repr(backend)


# ===================================================================
# K. Pipeline ↔ Backend integration
# ===================================================================


class TestPipelineBackendIntegration:
    """Integration tests verifying pipeline context carries backend info."""

    def test_context_holds_backend(self, default_config, methane):
        """Backend can be set on context and is accessible by stages."""
        ctx = PipelineContext(config=default_config, molecule=methane)
        backend = select_backend(default_config)
        ctx.backend = backend
        assert ctx.backend is backend
        assert ctx.backend.name == "pyscf"

    def test_stage_fails_without_backend(self, context_with_molecule):
        """Geometry opt stage requires both molecule and backend."""
        stage = GeometryOptimizationStage()
        result = stage.execute(context_with_molecule)
        assert result.status is StageStatus.FAILED
        assert "No QM backend" in result.message

    def test_input_stage_succeeds_without_backend(
        self, context_with_molecule
    ):
        """Input preparation doesn't need a backend."""
        stage = InputPreparationStage()
        result = stage.execute(context_with_molecule)
        assert result.status is StageStatus.COMPLETED

    def test_runner_with_single_stage_and_backend(
        self, default_config, methane
    ):
        """Pipeline runner + input stage succeeds with molecule set."""
        ctx = PipelineContext(config=default_config, molecule=methane)
        runner = PipelineRunner(stages=[InputPreparationStage()])
        result = runner.run(ctx)
        assert (
            result.stage_results["input_preparation"].status
            is StageStatus.COMPLETED
        )

    def test_backend_selection_from_config_variants(self):
        """select_backend produces the right type for each config variant."""
        configs = [
            (PoltypeConfig(), PySCFBackend),
            (PoltypeConfig(qm=QMConfig(use_gaus=True)), GaussianBackend),
            (
                PoltypeConfig(qm=QMConfig(dont_use_pyscf=True)),
                Psi4Backend,
            ),
        ]
        for cfg, expected_type in configs:
            backend = select_backend(cfg)
            assert isinstance(backend, expected_type), (
                f"Expected {expected_type.__name__} for config, "
                f"got {type(backend).__name__}"
            )

    def test_finalization_stage_without_backend(self, context):
        """Finalization doesn't require backend — always succeeds."""
        result = FinalizationStage().execute(context)
        assert result.status is StageStatus.COMPLETED

    def test_pipeline_context_artifact_flow(self, default_config, methane):
        """Artifacts set by one stage are visible to subsequent stages."""
        ctx = PipelineContext(config=default_config, molecule=methane)
        runner = PipelineRunner(
            stages=[
                _ArtifactStage("producer", "data", [1, 2, 3]),
                _CountingStage("consumer"),
            ]
        )
        result = runner.run(ctx)
        assert result.get_artifact("data") == [1, 2, 3]
        assert "producer" in result.stage_results
        assert "consumer" in result.stage_results
