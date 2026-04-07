"""
tests/unit/test_phase5.py – unit tests for Phase 5: checkpoint/resume & events.

Tests cover:
- PipelineEvent enum values.
- EventData construction and auto-timestamp.
- EventBus registration, emission ordering, and error-resilient dispatch.
- CheckpointData default construction.
- CheckpointManager save / load round-trip, should_skip_stage, and clear.
- PipelineRunner integration with EventBus (event emission + ordering).
- PipelineRunner integration with CheckpointManager (resume behaviour).
- CLI --resume and --no-checkpoint argument parsing.
- CLI _logging_hook dispatch.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional, Tuple
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from poltype.config.schema import PoltypeConfig
from poltype.molecule.molecule import Molecule
from poltype.pipeline.checkpoint import (
    CHECKPOINT_FILENAME,
    CheckpointData,
    CheckpointManager,
)
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.events import EventBus, EventData, PipelineEvent
from poltype.pipeline.factory import build_default_pipeline
from poltype.pipeline.runner import PipelineRunner
from poltype.pipeline.stage import Stage, StageResult, StageStatus
from poltype.qm.backend import (
    DMAResult,
    ESPGridResult,
    OptimizationResult,
    QMBackend,
    TorsionScanResult,
    WBOResult,
)


# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------


class PassStage(Stage):
    """Stage that always completes successfully."""

    def execute(self, context: PipelineContext) -> StageResult:
        return StageResult(
            status=StageStatus.COMPLETED,
            message=f"Stage '{self.name}' passed",
            artifacts={self.name: True},
        )


class FailStage(Stage):
    """Stage that always fails."""

    def execute(self, context: PipelineContext) -> StageResult:
        return StageResult(
            status=StageStatus.FAILED,
            message=f"Stage '{self.name}' failed on purpose",
        )


class MockBackend(QMBackend):
    """In-memory mock backend reused from Phase 4 tests."""

    @property
    def name(self) -> str:
        return "mock"

    def optimize_geometry(self, molecule, method="MP2", basis_set="6-31G*", **kw):
        return OptimizationResult(
            coordinates=molecule.coordinates, energy=-76.0, converged=True,
        )

    def compute_esp_grid(self, molecule, method="MP2", basis_set="aug-cc-pVTZ", **kw):
        return ESPGridResult(
            grid_coords=np.zeros((10, 3)), esp_values=np.zeros(10),
        )

    def torsion_scan(
        self, molecule, rotatable_bond, method="wB97X-D", basis_set="6-31G*",
        angles=None, **kw,
    ):
        scan = angles if angles is not None else np.arange(-180, 180, 15.0)
        return TorsionScanResult(angles=scan, energies=np.zeros_like(scan))

    def compute_dma(self, molecule, method="MP2", basis_set="6-311G**", **kw):
        from pathlib import Path
        fchk = Path(molecule.work_dir) / f"{molecule.name or 'mol'}_dma.fchk"
        fchk.write_text("mock fchk")
        return DMAResult(fchk_path=fchk, energy=-76.0)

    def compute_wbo_matrix(self, molecule, method="HF", basis_set="6-31G*", **kw):
        return WBOResult(wbo_matrix=np.eye(molecule.num_atoms))


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
def mock_backend():
    return MockBackend()


# ===================================================================
# A. PipelineEvent enum
# ===================================================================


class TestPipelineEvent:
    def test_has_expected_members(self):
        names = {e.name for e in PipelineEvent}
        assert "PIPELINE_STARTED" in names
        assert "PIPELINE_COMPLETED" in names
        assert "STAGE_STARTED" in names
        assert "STAGE_COMPLETED" in names
        assert "STAGE_SKIPPED" in names
        assert "STAGE_FAILED" in names

    def test_values_are_strings(self):
        for e in PipelineEvent:
            assert isinstance(e.value, str)

    def test_member_count(self):
        assert len(PipelineEvent) == 6


# ===================================================================
# B. EventData dataclass
# ===================================================================


class TestEventData:
    def test_default_construction(self):
        data = EventData(event=PipelineEvent.PIPELINE_STARTED)
        assert data.event is PipelineEvent.PIPELINE_STARTED
        assert data.stage_name == ""
        assert data.message == ""
        assert data.timestamp > 0
        assert data.extra == {}

    def test_custom_fields(self):
        data = EventData(
            event=PipelineEvent.STAGE_COMPLETED,
            stage_name="geometry_optimization",
            message="done",
            extra={"elapsed_seconds": 1.23},
        )
        assert data.stage_name == "geometry_optimization"
        assert data.extra["elapsed_seconds"] == 1.23

    def test_auto_timestamp_nonzero(self):
        d1 = EventData(event=PipelineEvent.PIPELINE_STARTED)
        d2 = EventData(event=PipelineEvent.PIPELINE_COMPLETED)
        assert d1.timestamp > 0
        assert d2.timestamp >= d1.timestamp

    def test_explicit_timestamp_preserved(self):
        data = EventData(
            event=PipelineEvent.PIPELINE_STARTED,
            timestamp=42.0,
        )
        assert data.timestamp == 42.0


# ===================================================================
# C. EventBus
# ===================================================================


class TestEventBus:
    def test_empty_bus_does_not_crash(self):
        bus = EventBus()
        bus.emit(EventData(event=PipelineEvent.PIPELINE_STARTED))

    def test_register_and_emit(self):
        received: List[EventData] = []
        bus = EventBus()
        bus.register(lambda d: received.append(d))
        data = EventData(event=PipelineEvent.STAGE_STARTED, stage_name="a")
        bus.emit(data)
        assert len(received) == 1
        assert received[0] is data

    def test_hooks_called_in_order(self):
        order: List[str] = []
        bus = EventBus()
        bus.register(lambda d: order.append("first"))
        bus.register(lambda d: order.append("second"))
        bus.register(lambda d: order.append("third"))
        bus.emit(EventData(event=PipelineEvent.PIPELINE_STARTED))
        assert order == ["first", "second", "third"]

    def test_multiple_hooks_receive_same_event(self):
        counts = [0, 0]

        def hook0(d):
            counts[0] += 1

        def hook1(d):
            counts[1] += 1

        bus = EventBus(hooks=[hook0, hook1])
        bus.emit(EventData(event=PipelineEvent.PIPELINE_STARTED))
        assert counts == [1, 1]

    def test_failing_hook_does_not_abort(self):
        received: List[str] = []

        def bad_hook(d):
            raise ValueError("boom")

        def good_hook(d):
            received.append("ok")

        bus = EventBus()
        bus.register(bad_hook)
        bus.register(good_hook)
        bus.emit(EventData(event=PipelineEvent.PIPELINE_STARTED))
        # good_hook should still be called despite bad_hook raising
        assert received == ["ok"]

    def test_hooks_property_returns_copy(self):
        bus = EventBus()

        def hook(d):
            pass

        bus.register(hook)
        hooks = bus.hooks
        hooks.clear()
        assert len(bus.hooks) == 1

    def test_repr(self):
        bus = EventBus()
        assert "EventBus" in repr(bus)
        assert "0" in repr(bus)

    def test_initial_hooks_in_constructor(self):
        hooks = [lambda d: None, lambda d: None]
        bus = EventBus(hooks=hooks)
        assert len(bus.hooks) == 2


# ===================================================================
# D. CheckpointData
# ===================================================================


class TestCheckpointData:
    def test_default_construction(self):
        data = CheckpointData()
        assert data.run_id == ""
        assert data.completed_stages == []
        assert data.stage_statuses == {}
        assert data.stage_messages == {}
        assert data.timestamp == 0.0

    def test_custom_construction(self):
        data = CheckpointData(
            run_id="abc123",
            completed_stages=["a", "b"],
            stage_statuses={"a": "completed", "b": "completed"},
        )
        assert data.run_id == "abc123"
        assert len(data.completed_stages) == 2

    def test_instances_are_independent(self):
        d1 = CheckpointData()
        d2 = CheckpointData()
        d1.completed_stages.append("x")
        assert d2.completed_stages == []


# ===================================================================
# E. CheckpointManager
# ===================================================================


class TestCheckpointManager:
    def test_checkpoint_path(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        assert mgr.checkpoint_path == tmp_path / CHECKPOINT_FILENAME

    def test_custom_filename(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path, filename="custom.json")
        assert mgr.checkpoint_path.name == "custom.json"

    def test_save_creates_file(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        result = StageResult(status=StageStatus.COMPLETED, message="ok")
        mgr.save("stage_a", result)
        assert mgr.checkpoint_path.exists()

    def test_save_records_completed_stage(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        result = StageResult(status=StageStatus.COMPLETED, message="ok")
        mgr.save("stage_a", result)
        assert "stage_a" in mgr.data.completed_stages
        assert mgr.data.stage_statuses["stage_a"] == "completed"

    def test_save_records_failed_stage(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        result = StageResult(status=StageStatus.FAILED, message="crash")
        mgr.save("stage_b", result)
        assert "stage_b" not in mgr.data.completed_stages
        assert mgr.data.stage_statuses["stage_b"] == "failed"

    def test_save_records_skipped_stage(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        result = StageResult(status=StageStatus.SKIPPED, message="skip")
        mgr.save("stage_c", result)
        assert "stage_c" not in mgr.data.completed_stages
        assert mgr.data.stage_statuses["stage_c"] == "skipped"

    def test_should_skip_stage(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        assert not mgr.should_skip_stage("stage_a")
        result = StageResult(status=StageStatus.COMPLETED, message="ok")
        mgr.save("stage_a", result)
        assert mgr.should_skip_stage("stage_a")
        assert not mgr.should_skip_stage("stage_b")

    def test_load_roundtrip(self, tmp_path):
        mgr1 = CheckpointManager(checkpoint_dir=tmp_path)
        mgr1.save("a", StageResult(status=StageStatus.COMPLETED, message="ok_a"))
        mgr1.save("b", StageResult(status=StageStatus.COMPLETED, message="ok_b"))
        mgr1.save("c", StageResult(status=StageStatus.FAILED, message="bad_c"))

        # Load into a fresh manager
        mgr2 = CheckpointManager(checkpoint_dir=tmp_path)
        data = mgr2.load()
        assert data.completed_stages == ["a", "b"]
        assert data.stage_statuses["c"] == "failed"
        assert data.stage_messages["a"] == "ok_a"
        assert data.stage_messages["c"] == "bad_c"

    def test_load_missing_file_raises(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        with pytest.raises(FileNotFoundError):
            mgr.load()

    def test_clear_removes_file(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        mgr.save("x", StageResult(status=StageStatus.COMPLETED))
        assert mgr.checkpoint_path.exists()
        mgr.clear()
        assert not mgr.checkpoint_path.exists()
        assert mgr.data.completed_stages == []

    def test_clear_no_file_does_not_raise(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        mgr.clear()  # should not raise

    def test_duplicate_save_does_not_duplicate_completed(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        result = StageResult(status=StageStatus.COMPLETED, message="ok")
        mgr.save("a", result)
        mgr.save("a", result)
        assert mgr.data.completed_stages.count("a") == 1

    def test_checkpoint_json_is_valid(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        mgr.save("a", StageResult(status=StageStatus.COMPLETED, message="ok"))
        raw = json.loads(mgr.checkpoint_path.read_text())
        assert "run_id" in raw
        assert "completed_stages" in raw
        assert "stage_statuses" in raw
        assert "stage_messages" in raw
        assert "timestamp" in raw

    def test_repr(self, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        r = repr(mgr)
        assert "CheckpointManager" in r
        assert "completed=0" in r

    def test_load_restores_run_id(self, tmp_path):
        mgr1 = CheckpointManager(checkpoint_dir=tmp_path)
        original_run_id = mgr1.data.run_id
        mgr1.save("a", StageResult(status=StageStatus.COMPLETED))

        mgr2 = CheckpointManager(checkpoint_dir=tmp_path)
        data = mgr2.load()
        assert data.run_id == original_run_id


# ===================================================================
# F. PipelineRunner + EventBus integration
# ===================================================================


class TestRunnerEventIntegration:
    def test_events_emitted_for_passing_stages(self, default_config):
        events: List[EventData] = []
        bus = EventBus(hooks=[lambda d: events.append(d)])

        runner = PipelineRunner(event_bus=bus)
        runner.add_stage(PassStage("a"))
        runner.add_stage(PassStage("b"))

        ctx = PipelineContext(config=default_config)
        runner.run(ctx)

        event_types = [e.event for e in events]
        assert event_types[0] is PipelineEvent.PIPELINE_STARTED
        assert PipelineEvent.STAGE_STARTED in event_types
        assert PipelineEvent.STAGE_COMPLETED in event_types
        assert event_types[-1] is PipelineEvent.PIPELINE_COMPLETED

    def test_stage_started_before_completed(self, default_config):
        events: List[EventData] = []
        bus = EventBus(hooks=[lambda d: events.append(d)])

        runner = PipelineRunner(event_bus=bus)
        runner.add_stage(PassStage("a"))
        ctx = PipelineContext(config=default_config)
        runner.run(ctx)

        stage_events = [e for e in events if e.stage_name == "a"]
        assert len(stage_events) == 2
        assert stage_events[0].event is PipelineEvent.STAGE_STARTED
        assert stage_events[1].event is PipelineEvent.STAGE_COMPLETED

    def test_failed_stage_emits_stage_failed(self, default_config):
        events: List[EventData] = []
        bus = EventBus(hooks=[lambda d: events.append(d)])

        runner = PipelineRunner(event_bus=bus)
        runner.add_stage(PassStage("a"))
        runner.add_stage(FailStage("b"))
        runner.add_stage(PassStage("c"))

        ctx = PipelineContext(config=default_config)
        runner.run(ctx)

        event_types = [e.event for e in events]
        assert PipelineEvent.STAGE_FAILED in event_types
        # Stage c should NOT have started or completed
        c_events = [e for e in events if e.stage_name == "c"]
        assert len(c_events) == 0

    def test_skipped_stage_emits_skipped(self, default_config):
        events: List[EventData] = []
        bus = EventBus(hooks=[lambda d: events.append(d)])

        class SkipStage(Stage):
            def execute(self, context):
                return StageResult(status=StageStatus.COMPLETED)

            def should_skip(self, context):
                return True

        runner = PipelineRunner(event_bus=bus)
        runner.add_stage(SkipStage("s"))
        ctx = PipelineContext(config=default_config)
        runner.run(ctx)

        skip_events = [
            e for e in events if e.event is PipelineEvent.STAGE_SKIPPED
        ]
        assert len(skip_events) == 1
        assert skip_events[0].stage_name == "s"

    def test_pipeline_completed_event_always_emitted(self, default_config):
        events: List[EventData] = []
        bus = EventBus(hooks=[lambda d: events.append(d)])

        runner = PipelineRunner(event_bus=bus)
        runner.add_stage(FailStage("x"))

        ctx = PipelineContext(config=default_config)
        runner.run(ctx)

        assert events[-1].event is PipelineEvent.PIPELINE_COMPLETED

    def test_event_bus_property(self):
        runner = PipelineRunner()
        assert runner.event_bus is None
        bus = EventBus()
        runner.event_bus = bus
        assert runner.event_bus is bus

    def test_no_event_bus_does_not_crash(self, default_config):
        runner = PipelineRunner()
        runner.add_stage(PassStage("a"))
        ctx = PipelineContext(config=default_config)
        result = runner.run(ctx)
        assert "a" in result.stage_results

    def test_elapsed_seconds_in_completed_event_extra(self, default_config):
        events: List[EventData] = []
        bus = EventBus(hooks=[lambda d: events.append(d)])

        runner = PipelineRunner(event_bus=bus)
        runner.add_stage(PassStage("a"))
        ctx = PipelineContext(config=default_config)
        runner.run(ctx)

        completed = [
            e for e in events if e.event is PipelineEvent.STAGE_COMPLETED
        ]
        assert len(completed) == 1
        assert "elapsed_seconds" in completed[0].extra


# ===================================================================
# G. PipelineRunner + CheckpointManager integration
# ===================================================================


class TestRunnerCheckpointIntegration:
    def test_checkpoint_saves_after_each_stage(self, default_config, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        runner = PipelineRunner(checkpoint_manager=mgr)
        runner.add_stage(PassStage("a"))
        runner.add_stage(PassStage("b"))

        ctx = PipelineContext(config=default_config)
        runner.run(ctx)

        assert mgr.should_skip_stage("a")
        assert mgr.should_skip_stage("b")

    def test_resume_skips_completed_stages(self, default_config, tmp_path):
        # First run: complete stages a and b, then fail on c
        mgr1 = CheckpointManager(checkpoint_dir=tmp_path)
        runner1 = PipelineRunner(checkpoint_manager=mgr1)
        runner1.add_stage(PassStage("a"))
        runner1.add_stage(PassStage("b"))
        runner1.add_stage(FailStage("c"))
        ctx1 = PipelineContext(config=default_config)
        runner1.run(ctx1)

        # Second run: resume should skip a and b
        mgr2 = CheckpointManager(checkpoint_dir=tmp_path)
        mgr2.load()

        events: List[EventData] = []
        bus = EventBus(hooks=[lambda d: events.append(d)])

        runner2 = PipelineRunner(event_bus=bus, checkpoint_manager=mgr2)
        runner2.add_stage(PassStage("a"))
        runner2.add_stage(PassStage("b"))
        runner2.add_stage(PassStage("c"))  # now passes
        ctx2 = PipelineContext(config=default_config)
        result = runner2.run(ctx2)

        # a and b should be skipped (checkpoint resume)
        assert result.stage_results["a"].status is StageStatus.SKIPPED
        assert "previous run" in result.stage_results["a"].message
        assert result.stage_results["b"].status is StageStatus.SKIPPED
        # c should complete normally
        assert result.stage_results["c"].status is StageStatus.COMPLETED

    def test_checkpoint_records_failed_stage(self, default_config, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        runner = PipelineRunner(checkpoint_manager=mgr)
        runner.add_stage(PassStage("a"))
        runner.add_stage(FailStage("b"))

        ctx = PipelineContext(config=default_config)
        runner.run(ctx)

        assert mgr.data.stage_statuses["a"] == "completed"
        assert mgr.data.stage_statuses["b"] == "failed"
        assert "b" not in mgr.data.completed_stages

    def test_no_checkpoint_manager_does_not_crash(self, default_config):
        runner = PipelineRunner()
        runner.add_stage(PassStage("a"))
        ctx = PipelineContext(config=default_config)
        result = runner.run(ctx)
        assert result.stage_results["a"].status is StageStatus.COMPLETED

    def test_checkpoint_manager_property(self, tmp_path):
        runner = PipelineRunner()
        assert runner.checkpoint_manager is None
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        runner.checkpoint_manager = mgr
        assert runner.checkpoint_manager is mgr

    def test_checkpoint_and_events_together(self, default_config, tmp_path):
        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        events: List[EventData] = []
        bus = EventBus(hooks=[lambda d: events.append(d)])

        runner = PipelineRunner(event_bus=bus, checkpoint_manager=mgr)
        runner.add_stage(PassStage("a"))
        runner.add_stage(PassStage("b"))

        ctx = PipelineContext(config=default_config)
        runner.run(ctx)

        # Both should be completed
        assert mgr.should_skip_stage("a")
        assert mgr.should_skip_stage("b")
        # Events should include both stage completed events
        completed_events = [
            e for e in events if e.event is PipelineEvent.STAGE_COMPLETED
        ]
        assert len(completed_events) == 2

    def test_skipped_stages_saved_in_checkpoint(self, default_config, tmp_path):
        class SkipStage(Stage):
            def execute(self, context):
                return StageResult(status=StageStatus.COMPLETED)

            def should_skip(self, context):
                return True

        mgr = CheckpointManager(checkpoint_dir=tmp_path)
        runner = PipelineRunner(checkpoint_manager=mgr)
        runner.add_stage(SkipStage("s"))
        ctx = PipelineContext(config=default_config)
        runner.run(ctx)

        assert mgr.data.stage_statuses["s"] == "skipped"


# ===================================================================
# H. Full pipeline with checkpoint resume
# ===================================================================


class TestFullPipelineWithCheckpoint:
    def test_full_pipeline_with_events_and_checkpoint(
        self, default_config, methane, mock_backend, tmp_path,
    ):
        events: List[EventData] = []
        bus = EventBus(hooks=[lambda d: events.append(d)])
        mgr = CheckpointManager(checkpoint_dir=tmp_path)

        runner = build_default_pipeline()
        runner.event_bus = bus
        runner.checkpoint_manager = mgr

        ctx = PipelineContext(
            config=default_config,
            molecule=methane,
            backend=mock_backend,
        )
        result = runner.run(ctx)

        # All 10 stages complete
        assert len(result.stage_results) == 10
        for sr in result.stage_results.values():
            assert sr.status is StageStatus.COMPLETED

        # Checkpoint has all 10 stages
        assert len(mgr.data.completed_stages) == 10

        # Events include pipeline start/end + 10 started + 10 completed = 22
        assert events[0].event is PipelineEvent.PIPELINE_STARTED
        assert events[-1].event is PipelineEvent.PIPELINE_COMPLETED
        started = [e for e in events if e.event is PipelineEvent.STAGE_STARTED]
        completed = [e for e in events if e.event is PipelineEvent.STAGE_COMPLETED]
        assert len(started) == 10
        assert len(completed) == 10

    def test_resume_full_pipeline(
        self, default_config, methane, mock_backend, tmp_path,
    ):
        # First run
        mgr1 = CheckpointManager(checkpoint_dir=tmp_path)
        runner1 = build_default_pipeline()
        runner1.checkpoint_manager = mgr1
        ctx1 = PipelineContext(
            config=default_config,
            molecule=methane,
            backend=mock_backend,
        )
        runner1.run(ctx1)

        # Second run: all stages should be skipped on resume
        mgr2 = CheckpointManager(checkpoint_dir=tmp_path)
        mgr2.load()

        runner2 = build_default_pipeline()
        runner2.checkpoint_manager = mgr2

        ctx2 = PipelineContext(
            config=default_config,
            molecule=methane,
            backend=mock_backend,
        )
        result2 = runner2.run(ctx2)

        for sr in result2.stage_results.values():
            assert sr.status is StageStatus.SKIPPED


# ===================================================================
# I. CLI argument parsing
# ===================================================================


class TestCLIPhase5Args:
    def test_resume_flag(self):
        from poltype.__main__ import _build_parser

        parser = _build_parser()
        args = parser.parse_args(["--resume"])
        assert args.resume is True

    def test_no_checkpoint_flag(self):
        from poltype.__main__ import _build_parser

        parser = _build_parser()
        args = parser.parse_args(["--no-checkpoint"])
        assert args.no_checkpoint is True

    def test_resume_default_false(self):
        from poltype.__main__ import _build_parser

        parser = _build_parser()
        args = parser.parse_args([])
        assert args.resume is False

    def test_no_checkpoint_default_false(self):
        from poltype.__main__ import _build_parser

        parser = _build_parser()
        args = parser.parse_args([])
        assert args.no_checkpoint is False

    def test_resume_and_no_checkpoint_together(self):
        from poltype.__main__ import _build_parser

        parser = _build_parser()
        args = parser.parse_args(["--resume", "--no-checkpoint"])
        assert args.resume is True
        assert args.no_checkpoint is True


# ===================================================================
# J. CLI _logging_hook
# ===================================================================


class TestLoggingHook:
    def test_logging_hook_handles_all_events(self):
        from poltype.__main__ import _logging_hook

        # Should not raise for any event type
        for event in PipelineEvent:
            data = EventData(
                event=event,
                stage_name="test_stage",
                message="test message",
                extra={"elapsed_seconds": 1.0},
            )
            _logging_hook(data)  # should not raise

    def test_logging_hook_is_callable(self):
        from poltype.__main__ import _logging_hook

        assert callable(_logging_hook)


# ===================================================================
# K. Pipeline __init__ exports
# ===================================================================


class TestPipelineExports:
    def test_event_types_importable(self):
        from poltype.pipeline import EventBus, EventData, PipelineEvent

        assert EventBus is not None
        assert EventData is not None
        assert PipelineEvent is not None

    def test_checkpoint_types_importable(self):
        from poltype.pipeline import CheckpointData, CheckpointManager

        assert CheckpointData is not None
        assert CheckpointManager is not None
