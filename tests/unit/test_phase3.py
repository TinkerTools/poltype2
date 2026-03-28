"""
tests/unit/test_phase3.py – unit tests for Phase 3: CLI, logging,
pipeline factory, and error hierarchy.

Tests cover:
- Custom exception hierarchy (PoltypeError, ConfigError, StageError,
  BackendError).
- Logging configuration (setup_logging).
- Pipeline factory (build_default_pipeline).
- CLI argument parsing and main() entry point.
"""

from __future__ import annotations

import logging
from pathlib import Path
from unittest.mock import patch

import pytest

from poltype.config.schema import PoltypeConfig
from poltype.errors import BackendError, ConfigError, PoltypeError, StageError
from poltype.logging_config import LOGGER_NAME, setup_logging
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.factory import build_default_pipeline
from poltype.pipeline.stage import StageStatus


# ---------------------------------------------------------------------------
# Error hierarchy
# ---------------------------------------------------------------------------


class TestPoltypeError:
    def test_is_exception(self):
        assert issubclass(PoltypeError, Exception)

    def test_instantiation(self):
        err = PoltypeError("something went wrong")
        assert str(err) == "something went wrong"


class TestConfigError:
    def test_is_poltype_error(self):
        assert issubclass(ConfigError, PoltypeError)

    def test_caught_by_base(self):
        with pytest.raises(PoltypeError):
            raise ConfigError("bad config")


class TestStageError:
    def test_is_poltype_error(self):
        assert issubclass(StageError, PoltypeError)

    def test_stage_name_attribute(self):
        err = StageError("failed", stage_name="geometry_optimization")
        assert err.stage_name == "geometry_optimization"
        assert str(err) == "failed"

    def test_default_stage_name(self):
        err = StageError("failed")
        assert err.stage_name == ""


class TestBackendError:
    def test_is_poltype_error(self):
        assert issubclass(BackendError, PoltypeError)

    def test_backend_name_attribute(self):
        err = BackendError("psi4 crashed", backend_name="psi4")
        assert err.backend_name == "psi4"
        assert str(err) == "psi4 crashed"

    def test_default_backend_name(self):
        err = BackendError("unknown error")
        assert err.backend_name == ""


# ---------------------------------------------------------------------------
# Logging configuration
# ---------------------------------------------------------------------------


class TestSetupLogging:
    def setup_method(self):
        """Remove all handlers from the poltype logger before each test."""
        logger = logging.getLogger(LOGGER_NAME)
        logger.handlers.clear()

    def test_returns_logger(self):
        result = setup_logging()
        assert isinstance(result, logging.Logger)
        assert result.name == LOGGER_NAME

    def test_default_level_is_info(self):
        logger = setup_logging()
        assert logger.level == logging.INFO

    def test_verbose_sets_debug(self):
        logger = setup_logging(verbose=True)
        assert logger.level == logging.DEBUG

    def test_quiet_sets_warning(self):
        logger = setup_logging(quiet=True)
        assert logger.level == logging.WARNING

    def test_quiet_overrides_verbose(self):
        logger = setup_logging(verbose=True, quiet=True)
        assert logger.level == logging.WARNING

    def test_adds_stderr_handler(self):
        logger = setup_logging()
        assert len(logger.handlers) >= 1
        assert isinstance(logger.handlers[0], logging.StreamHandler)

    def test_log_file_handler(self, tmp_path):
        log_file = str(tmp_path / "test.log")
        logger = setup_logging(log_file=log_file)
        assert len(logger.handlers) == 2
        assert isinstance(logger.handlers[1], logging.FileHandler)

    def test_no_duplicate_handlers_on_repeated_calls(self):
        setup_logging()
        setup_logging()
        logger = logging.getLogger(LOGGER_NAME)
        # Second call should update levels, not add new handlers
        assert len(logger.handlers) == 1


# ---------------------------------------------------------------------------
# Pipeline factory
# ---------------------------------------------------------------------------


class TestBuildDefaultPipeline:
    def test_returns_pipeline_runner(self):
        from poltype.pipeline.runner import PipelineRunner

        runner = build_default_pipeline()
        assert isinstance(runner, PipelineRunner)

    def test_has_six_stages(self):
        runner = build_default_pipeline()
        assert repr(runner) == "PipelineRunner(stages=6)"

    def test_stage_order(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert names == [
            "input_preparation",
            "geometry_optimization",
            "esp_fitting",
            "multipole",
            "torsion_fitting",
            "finalization",
        ]

    def test_runs_with_molecule_context(self):
        """The default pipeline should run to completion with a molecule."""
        from poltype.molecule.molecule import Molecule

        runner = build_default_pipeline()
        mol = Molecule.from_smiles("C", name="methane")
        cfg = PoltypeConfig()
        ctx = PipelineContext(config=cfg, molecule=mol)
        result = runner.run(ctx)
        # All 6 stages should have recorded a result
        assert len(result.stage_results) == 6
        for sr in result.stage_results.values():
            assert sr.status is StageStatus.COMPLETED


# ---------------------------------------------------------------------------
# Pipeline factory re-export
# ---------------------------------------------------------------------------


class TestPipelinePackageReexport:
    def test_build_default_pipeline_importable(self):
        from poltype.pipeline import build_default_pipeline as bpf

        assert callable(bpf)


# ---------------------------------------------------------------------------
# Top-level error re-exports
# ---------------------------------------------------------------------------


class TestTopLevelReexports:
    def test_poltype_error_importable(self):
        from poltype import PoltypeError

        assert issubclass(PoltypeError, Exception)

    def test_config_error_importable(self):
        from poltype import ConfigError

        assert issubclass(ConfigError, PoltypeError)

    def test_stage_error_importable(self):
        from poltype import StageError

        assert issubclass(StageError, PoltypeError)

    def test_backend_error_importable(self):
        from poltype import BackendError

        assert issubclass(BackendError, PoltypeError)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


class TestCLI:
    def test_version_flag(self, capsys):
        from poltype.__main__ import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0

    def test_help_flag(self, capsys):
        from poltype.__main__ import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--help"])
        assert exc_info.value.code == 0

    def test_missing_config_returns_2(self, tmp_path):
        from poltype.__main__ import main

        code = main(["--config", str(tmp_path / "nonexistent.ini")])
        assert code == 2

    def test_successful_run(self, tmp_path):
        """End-to-end: config → pipeline → success exit code."""
        from poltype.__main__ import main

        ini = tmp_path / "poltype.ini"
        ini.write_text("structure = dummy.sdf\n")

        # Create a minimal SDF file for the molecule loader
        sdf = tmp_path / "dummy.sdf"
        from rdkit import Chem

        mol = Chem.MolFromSmiles("C")
        mol = Chem.AddHs(mol)
        writer = Chem.SDWriter(str(sdf))
        writer.write(mol)
        writer.close()

        code = main([
            "--config", str(ini),
            "--work-dir", str(tmp_path),
            "--quiet",
        ])
        assert code == 0

    def test_verbose_flag_accepted(self, tmp_path):
        from poltype.__main__ import main

        ini = tmp_path / "poltype.ini"
        ini.write_text("structure = dummy.sdf\n")

        sdf = tmp_path / "dummy.sdf"
        from rdkit import Chem

        mol = Chem.MolFromSmiles("C")
        mol = Chem.AddHs(mol)
        writer = Chem.SDWriter(str(sdf))
        writer.write(mol)
        writer.close()

        code = main([
            "--config", str(ini),
            "--work-dir", str(tmp_path),
            "--verbose",
        ])
        assert code == 0


# ---------------------------------------------------------------------------
# Runner logging integration
# ---------------------------------------------------------------------------


class TestRunnerLogging:
    def setup_method(self):
        """Clear poltype logger handlers."""
        logger = logging.getLogger(LOGGER_NAME)
        logger.handlers.clear()

    def test_runner_emits_stage_logs(self, caplog):
        """The runner should emit INFO logs for each stage."""
        from poltype.molecule.molecule import Molecule
        from poltype.pipeline.runner import PipelineRunner
        from poltype.pipeline.stages.input_prep import InputPreparationStage

        setup_logging()
        mol = Molecule.from_smiles("C", name="methane")
        cfg = PoltypeConfig()
        ctx = PipelineContext(config=cfg, molecule=mol)

        runner = PipelineRunner(stages=[InputPreparationStage()])
        with caplog.at_level(logging.INFO, logger="poltype.pipeline.runner"):
            runner.run(ctx)

        messages = [r.message for r in caplog.records]
        assert any("starting" in m for m in messages)
        assert any("completed" in m for m in messages)

    def test_runner_logs_skip(self, caplog):
        """The runner should log when a stage is skipped."""
        from poltype.pipeline.runner import PipelineRunner
        from poltype.pipeline.stages.geometry_opt import (
            GeometryOptimizationStage,
        )

        setup_logging()
        cfg = PoltypeConfig()
        cfg.database_match_only = True
        ctx = PipelineContext(config=cfg)

        runner = PipelineRunner(stages=[GeometryOptimizationStage()])
        with caplog.at_level(logging.INFO, logger="poltype.pipeline.runner"):
            runner.run(ctx)

        messages = [r.message for r in caplog.records]
        assert any("skipped" in m for m in messages)

    def test_runner_logs_failure(self, caplog):
        """The runner should log ERROR when a stage fails."""
        from poltype.pipeline.context import PipelineContext
        from poltype.pipeline.runner import PipelineRunner
        from poltype.pipeline.stage import Stage, StageResult, StageStatus

        class _FailStage(Stage):
            def execute(self, context):
                return StageResult(
                    status=StageStatus.FAILED, message="boom"
                )

        setup_logging()
        cfg = PoltypeConfig()
        ctx = PipelineContext(config=cfg)

        runner = PipelineRunner(stages=[_FailStage("fail_stage")])
        with caplog.at_level(logging.ERROR, logger="poltype.pipeline.runner"):
            runner.run(ctx)

        messages = [r.message for r in caplog.records]
        assert any("failed" in m for m in messages)
