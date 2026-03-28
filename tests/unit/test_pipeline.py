"""
tests/unit/test_pipeline.py – unit tests for the pipeline framework.

Tests cover:
- StageStatus enum values.
- StageResult dataclass construction and defaults.
- Stage ABC enforcement and concrete subclass behaviour.
- PipelineContext construction, artifact helpers, and molecule updates.
- PipelineRunner ordering, failure halting, skip handling, and artifact
  merging.
- Concrete stage names, execute() stubs, and should_skip() logic.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest

from poltype.config.schema import PoltypeConfig
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
        # Should never be reached
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
# A. StageStatus tests
# ===================================================================


class TestStageStatus:
    def test_all_enum_values_exist(self):
        assert StageStatus.PENDING.value == "pending"
        assert StageStatus.RUNNING.value == "running"
        assert StageStatus.COMPLETED.value == "completed"
        assert StageStatus.SKIPPED.value == "skipped"
        assert StageStatus.FAILED.value == "failed"

    def test_enum_member_count(self):
        assert len(StageStatus) == 5


# ===================================================================
# B. StageResult tests
# ===================================================================


class TestStageResult:
    def test_instantiate_with_status_only(self):
        r = StageResult(status=StageStatus.COMPLETED)
        assert r.status is StageStatus.COMPLETED
        assert r.message == ""
        assert r.elapsed_seconds == 0.0
        assert r.artifacts == {}

    def test_instantiate_with_all_fields(self):
        r = StageResult(
            status=StageStatus.FAILED,
            message="boom",
            elapsed_seconds=1.23,
            artifacts={"key": "val"},
        )
        assert r.status is StageStatus.FAILED
        assert r.message == "boom"
        assert r.elapsed_seconds == 1.23
        assert r.artifacts == {"key": "val"}

    def test_artifacts_default_is_independent(self):
        """Each instance should get its own empty dict."""
        r1 = StageResult(status=StageStatus.COMPLETED)
        r2 = StageResult(status=StageStatus.COMPLETED)
        r1.artifacts["a"] = 1
        assert "a" not in r2.artifacts


# ===================================================================
# C. Stage ABC tests
# ===================================================================


class TestStageABC:
    def test_cannot_instantiate_directly(self):
        with pytest.raises(TypeError):
            Stage("some-stage")  # type: ignore[abstract]

    def test_concrete_subclass_works(self):
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
# D. PipelineContext tests
# ===================================================================


class TestPipelineContext:
    def test_instantiate_with_config_only(self, default_config):
        ctx = PipelineContext(config=default_config)
        assert ctx.config is default_config
        assert ctx.molecule is None
        assert ctx.backend is None
        assert ctx.work_dir == Path.cwd()

    def test_molecule_starts_none_when_not_provided(self, default_config):
        ctx = PipelineContext(config=default_config)
        assert ctx.molecule is None

    def test_work_dir_defaults_to_cwd(self, default_config):
        ctx = PipelineContext(config=default_config)
        assert ctx.work_dir == Path.cwd()

    def test_work_dir_custom(self, default_config, tmp_path):
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
        assert context.molecule.name == "methane"

    def test_stage_results_starts_empty(self, context):
        assert context.stage_results == {}

    def test_molecule_setter(self, context, methane):
        context.molecule = methane
        assert context.molecule is methane

    def test_repr(self, context):
        r = repr(context)
        assert "PipelineContext" in r
        assert "PoltypeConfig" in r

    def test_repr_with_molecule(self, context_with_molecule):
        r = repr(context_with_molecule)
        assert "methane" in r


# ===================================================================
# E. PipelineRunner tests
# ===================================================================


class TestPipelineRunner:
    def test_instantiate_empty(self):
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

    def test_run_empty_returns_context_unchanged(self, context):
        runner = PipelineRunner()
        result_ctx = runner.run(context)
        assert result_ctx is context
        assert result_ctx.stage_results == {}
        assert result_ctx.artifacts == {}

    def test_run_executes_stages_in_order(self, context):
        runner = PipelineRunner(
            stages=[
                _CountingStage("first"),
                _CountingStage("second"),
                _CountingStage("third"),
            ]
        )
        result_ctx = runner.run(context)
        assert result_ctx.get_artifact("order") == [
            "first",
            "second",
            "third",
        ]

    def test_run_records_stage_results(self, context):
        runner = PipelineRunner(
            stages=[_CountingStage("a"), _CountingStage("b")]
        )
        result_ctx = runner.run(context)
        assert "a" in result_ctx.stage_results
        assert "b" in result_ctx.stage_results
        assert result_ctx.stage_results["a"].status is StageStatus.COMPLETED
        assert result_ctx.stage_results["b"].status is StageStatus.COMPLETED

    def test_run_stops_after_failed_stage(self, context):
        runner = PipelineRunner(
            stages=[
                _CountingStage("before"),
                _FailingStage("boom"),
                _CountingStage("after"),
            ]
        )
        result_ctx = runner.run(context)
        assert "before" in result_ctx.stage_results
        assert "boom" in result_ctx.stage_results
        assert "after" not in result_ctx.stage_results
        assert (
            result_ctx.stage_results["boom"].status is StageStatus.FAILED
        )

    def test_run_skips_stages(self, context):
        runner = PipelineRunner(
            stages=[
                _CountingStage("before"),
                _AlwaysSkipStage("skipped"),
                _CountingStage("after"),
            ]
        )
        result_ctx = runner.run(context)
        assert (
            result_ctx.stage_results["skipped"].status is StageStatus.SKIPPED
        )
        # Skipped stage should not have put its artifact in context
        assert "should_not_appear" not in result_ctx.artifacts
        # Other stages still ran
        assert result_ctx.get_artifact("order") == ["before", "after"]

    def test_run_merges_artifacts(self, context):
        runner = PipelineRunner(
            stages=[
                _ArtifactStage("s1", "key1", "val1"),
                _ArtifactStage("s2", "key2", "val2"),
            ]
        )
        result_ctx = runner.run(context)
        assert result_ctx.get_artifact("key1") == "val1"
        assert result_ctx.get_artifact("key2") == "val2"

    def test_run_elapsed_seconds_populated(self, context):
        runner = PipelineRunner(stages=[_CountingStage("timed")])
        result_ctx = runner.run(context)
        assert result_ctx.stage_results["timed"].elapsed_seconds >= 0.0

    def test_instantiate_with_stages_list(self):
        stages = [_CountingStage("a"), _CountingStage("b")]
        runner = PipelineRunner(stages=stages)
        assert repr(runner) == "PipelineRunner(stages=2)"


# ===================================================================
# F. Concrete stage tests
# ===================================================================


class TestInputPreparationStage:
    def test_name(self):
        assert InputPreparationStage().name == "input_preparation"

    def test_execute_with_existing_molecule(self, context_with_molecule):
        stage = InputPreparationStage()
        result = stage.execute(context_with_molecule)
        assert result.status is StageStatus.COMPLETED
        assert "already present" in result.message

    def test_execute_no_structure_fails(self, context):
        stage = InputPreparationStage()
        result = stage.execute(context)
        assert result.status is StageStatus.FAILED
        assert "No structure" in result.message

    def test_should_skip_when_check_input_only(self, default_config):
        default_config.check_input_only = True
        ctx = PipelineContext(config=default_config)
        stage = InputPreparationStage()
        assert stage.should_skip(ctx) is True

    def test_should_not_skip_by_default(self, context):
        stage = InputPreparationStage()
        assert stage.should_skip(context) is False


class TestGeometryOptimizationStage:
    def test_name(self):
        assert GeometryOptimizationStage().name == "geometry_optimization"

    def test_execute_fails_when_no_molecule(self, context):
        stage = GeometryOptimizationStage()
        result = stage.execute(context)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_execute_fails_when_no_backend(self, context_with_molecule):
        stage = GeometryOptimizationStage()
        result = stage.execute(context_with_molecule)
        assert result.status is StageStatus.FAILED
        assert "No QM backend" in result.message

    def test_should_skip_when_database_match_only(self, default_config):
        default_config.database_match_only = True
        ctx = PipelineContext(config=default_config)
        stage = GeometryOptimizationStage()
        assert stage.should_skip(ctx) is True

    def test_should_not_skip_by_default(self, context):
        stage = GeometryOptimizationStage()
        assert stage.should_skip(context) is False


class TestESPFittingStage:
    def test_name(self):
        assert ESPFittingStage().name == "esp_fitting"

    def test_execute_fails_when_no_molecule(self, context):
        stage = ESPFittingStage()
        result = stage.execute(context)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_execute_fails_when_no_backend(self, context_with_molecule):
        stage = ESPFittingStage()
        result = stage.execute(context_with_molecule)
        assert result.status is StageStatus.FAILED
        assert "No QM backend" in result.message

    def test_should_skip_when_database_match_only(self, default_config):
        default_config.database_match_only = True
        ctx = PipelineContext(config=default_config)
        stage = ESPFittingStage()
        assert stage.should_skip(ctx) is True

    def test_should_not_skip_by_default(self, context):
        stage = ESPFittingStage()
        assert stage.should_skip(context) is False


class TestMultipoleStage:
    def test_name(self):
        assert MultipoleStage().name == "multipole"

    def test_execute_fails_when_no_molecule(self, context):
        stage = MultipoleStage()
        result = stage.execute(context)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_execute_completes_without_esp_data(self, context_with_molecule):
        stage = MultipoleStage()
        result = stage.execute(context_with_molecule)
        assert result.status is StageStatus.COMPLETED
        assert "No ESP data" in result.message

    def test_should_skip_when_database_match_only(self, default_config):
        default_config.database_match_only = True
        ctx = PipelineContext(config=default_config)
        stage = MultipoleStage()
        assert stage.should_skip(ctx) is True

    def test_should_not_skip_by_default(self, context):
        stage = MultipoleStage()
        assert stage.should_skip(context) is False


class TestTorsionFittingStage:
    def test_name(self):
        assert TorsionFittingStage().name == "torsion_fitting"

    def test_execute_fails_when_no_molecule(self, context):
        stage = TorsionFittingStage()
        result = stage.execute(context)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_execute_fails_when_no_backend(self, context_with_molecule):
        stage = TorsionFittingStage()
        result = stage.execute(context_with_molecule)
        assert result.status is StageStatus.FAILED
        assert "No QM backend" in result.message

    def test_should_skip_when_do_torsion_fit_false(self, default_config):
        default_config.do_torsion_fit = False
        ctx = PipelineContext(config=default_config)
        stage = TorsionFittingStage()
        assert stage.should_skip(ctx) is True

    def test_should_skip_when_database_match_only(self, default_config):
        default_config.database_match_only = True
        ctx = PipelineContext(config=default_config)
        stage = TorsionFittingStage()
        assert stage.should_skip(ctx) is True

    def test_should_not_skip_by_default(self, context):
        stage = TorsionFittingStage()
        assert stage.should_skip(context) is False


class TestFinalizationStage:
    def test_name(self):
        assert FinalizationStage().name == "finalization"

    def test_execute_returns_completed(self, context):
        stage = FinalizationStage()
        result = stage.execute(context)
        assert result.status is StageStatus.COMPLETED
        assert "final_summary" in result.artifacts

    def test_execute_collects_molecule_info(self, context_with_molecule):
        stage = FinalizationStage()
        result = stage.execute(context_with_molecule)
        assert result.status is StageStatus.COMPLETED
        summary = result.artifacts["final_summary"]
        assert summary["has_molecule"] is True
        assert summary["molecule_name"] == "methane"
