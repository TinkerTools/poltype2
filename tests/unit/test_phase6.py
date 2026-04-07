"""
tests/unit/test_phase6.py – unit tests for Phase 6: output writers & dry-run.

Tests cover:
- OutputWriter ABC interface enforcement.
- XYZFileWriter: writes valid Tinker XYZ with correct atom count, coords,
  neighbours.
- KeyFileWriter: writes AMOEBA key with parameter reference, atom defs,
  artifact annotations.
- RunSummary construction from pipeline context.
- SummaryWriter: produces human-readable summary file.
- FinalizationStage with writers: writes output files and records paths.
- FinalizationStage without writers: backwards-compatible summary.
- Dry-run mode: QM-heavy stages skip, input_prep and finalization run.
- Dry-run config field and loader integration.
- CLI --dry-run flag parsing.
- Full pipeline integration with dry-run.
"""

from __future__ import annotations

import textwrap
from pathlib import Path
from typing import List, Optional, Tuple
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from poltype.config.schema import PoltypeConfig, QMConfig
from poltype.molecule.molecule import Molecule
from poltype.output.key_writer import KeyFileWriter
from poltype.output.params import MultipoleParam, PolarizeParam, TorsionFold, TorsionParam
from poltype.output.summary import RunSummary, SummaryWriter
from poltype.output.writer import OutputWriter
from poltype.output.xyz_writer import XYZFileWriter
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.factory import build_default_pipeline
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


# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------


class MockBackend(QMBackend):
    """In-memory mock backend reused from Phase 4 tests."""

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
    """Create a simple ethanol molecule for testing."""
    return Molecule.from_smiles("CCO", name="ethanol")


def _make_context(
    tmp_path: Path,
    molecule: Optional[Molecule] = None,
    backend: Optional[QMBackend] = None,
    dry_run: bool = False,
) -> PipelineContext:
    """Create a pipeline context for testing."""
    config = PoltypeConfig(dry_run=dry_run)
    mol = molecule or _make_ethanol()
    return PipelineContext(
        config=config,
        molecule=mol,
        backend=backend or MockBackend(),
        work_dir=tmp_path,
    )


# ===================================================================
# OutputWriter ABC tests
# ===================================================================


class TestOutputWriterABC:
    """Tests for the OutputWriter abstract base class."""

    def test_cannot_instantiate_abc(self):
        """OutputWriter itself cannot be instantiated."""
        with pytest.raises(TypeError):
            OutputWriter()

    def test_concrete_must_implement_name_and_write(self):
        """A subclass missing name or write raises TypeError."""

        class Incomplete(OutputWriter):
            pass

        with pytest.raises(TypeError):
            Incomplete()

    def test_concrete_with_methods_can_instantiate(self):
        """A fully implemented subclass can be instantiated."""

        class DummyWriter(OutputWriter):
            @property
            def name(self):
                return "dummy"

            def write(self, context):
                return Path("/tmp/dummy.txt")

        w = DummyWriter()
        assert w.name == "dummy"

    def test_output_dir_property(self):
        """output_dir returns the configured directory."""

        class DummyWriter(OutputWriter):
            @property
            def name(self):
                return "dummy"

            def write(self, context):
                return Path("/tmp/dummy.txt")

        w = DummyWriter(output_dir=Path("/custom"))
        assert w.output_dir == Path("/custom")

    def test_output_dir_default_none(self):
        """output_dir defaults to None."""

        class DummyWriter(OutputWriter):
            @property
            def name(self):
                return "dummy"

            def write(self, context):
                return Path("/tmp/dummy.txt")

        w = DummyWriter()
        assert w.output_dir is None

    def test_repr(self):
        """repr includes the class name and output_dir."""

        class DummyWriter(OutputWriter):
            @property
            def name(self):
                return "dummy"

            def write(self, context):
                return Path("/tmp/dummy.txt")

        w = DummyWriter(output_dir=Path("/custom"))
        assert "DummyWriter" in repr(w)
        assert "/custom" in repr(w)

    def test_resolve_dir_falls_back_to_context(self, tmp_path):
        """_resolve_dir returns context.work_dir when output_dir is None."""

        class DummyWriter(OutputWriter):
            @property
            def name(self):
                return "dummy"

            def write(self, context):
                return self._resolve_dir(context)

        ctx = _make_context(tmp_path)
        w = DummyWriter()
        result = w.write(ctx)
        assert result == tmp_path

    def test_resolve_dir_uses_configured(self, tmp_path):
        """_resolve_dir uses configured output_dir when set."""
        custom = tmp_path / "custom_out"

        class DummyWriter(OutputWriter):
            @property
            def name(self):
                return "dummy"

            def write(self, context):
                return self._resolve_dir(context)

        ctx = _make_context(tmp_path)
        w = DummyWriter(output_dir=custom)
        result = w.write(ctx)
        assert result == custom
        assert custom.exists()


# ===================================================================
# XYZFileWriter tests
# ===================================================================


class TestXYZFileWriter:
    """Tests for the XYZFileWriter."""

    def test_name(self):
        w = XYZFileWriter()
        assert w.name == "xyz_writer"

    def test_write_creates_file(self, tmp_path):
        ctx = _make_context(tmp_path)
        w = XYZFileWriter()
        path = w.write(ctx)
        assert path.exists()
        assert path.name == "final.xyz"

    def test_write_custom_filename(self, tmp_path):
        ctx = _make_context(tmp_path)
        w = XYZFileWriter(filename="output.xyz")
        path = w.write(ctx)
        assert path.name == "output.xyz"
        assert path.exists()

    def test_write_custom_output_dir(self, tmp_path):
        ctx = _make_context(tmp_path)
        custom_dir = tmp_path / "xyz_output"
        w = XYZFileWriter(output_dir=custom_dir)
        path = w.write(ctx)
        assert path.parent == custom_dir

    def test_first_line_has_atom_count(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        w = XYZFileWriter()
        path = w.write(ctx)
        lines = path.read_text().strip().split("\n")
        first_line = lines[0].strip()
        # First line: "N  name"
        parts = first_line.split()
        assert int(parts[0]) == mol.num_atoms

    def test_first_line_has_molecule_name(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        w = XYZFileWriter()
        path = w.write(ctx)
        lines = path.read_text().strip().split("\n")
        assert "ethanol" in lines[0]

    def test_correct_number_of_atom_lines(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        w = XYZFileWriter()
        path = w.write(ctx)
        lines = path.read_text().strip().split("\n")
        # First line is header, rest are atom lines
        assert len(lines) == mol.num_atoms + 1

    def test_atom_lines_have_coordinates(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        w = XYZFileWriter()
        path = w.write(ctx)
        lines = path.read_text().strip().split("\n")
        # Check an atom line has at least 6 fields
        atom_line = lines[1].split()
        assert len(atom_line) >= 6  # idx, sym, x, y, z, type, [neighbours]

    def test_atom_indices_are_one_based(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        w = XYZFileWriter()
        path = w.write(ctx)
        lines = path.read_text().strip().split("\n")
        for i, line in enumerate(lines[1:], start=1):
            idx = int(line.split()[0])
            assert idx == i

    def test_raises_on_no_molecule(self, tmp_path):
        config = PoltypeConfig()
        ctx = PipelineContext(config=config, work_dir=tmp_path)
        w = XYZFileWriter()
        with pytest.raises(ValueError, match="No molecule"):
            w.write(ctx)

    def test_coordinates_match_molecule(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        w = XYZFileWriter()
        path = w.write(ctx)
        lines = path.read_text().strip().split("\n")
        coords = mol.coordinates
        for i, line in enumerate(lines[1:]):
            parts = line.split()
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            assert abs(x - coords[i][0]) < 1e-5
            assert abs(y - coords[i][1]) < 1e-5
            assert abs(z - coords[i][2]) < 1e-5


# ===================================================================
# KeyFileWriter tests
# ===================================================================


class TestKeyFileWriter:
    """Tests for the KeyFileWriter."""

    def test_name(self):
        w = KeyFileWriter()
        assert w.name == "key_writer"

    def test_write_creates_file(self, tmp_path):
        ctx = _make_context(tmp_path)
        w = KeyFileWriter()
        path = w.write(ctx)
        assert path.exists()
        assert path.name == "final.key"

    def test_write_custom_filename(self, tmp_path):
        ctx = _make_context(tmp_path)
        w = KeyFileWriter(filename="custom.key")
        path = w.write(ctx)
        assert path.name == "custom.key"

    def test_contains_parameters_line(self, tmp_path):
        ctx = _make_context(tmp_path)
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert "parameters" in content

    def test_contains_atom_definitions(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        # Should have 'atom' lines for each atom
        atom_lines = [l for l in content.split("\n") if l.startswith("atom")]
        assert len(atom_lines) == mol.num_atoms

    def test_contains_molecule_name_in_header(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert "ethanol" in content

    def test_contains_formula(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert mol.formula in content

    def test_opt_energy_annotation_when_present(self, tmp_path):
        ctx = _make_context(tmp_path)
        ctx.set_artifact(
            "opt_result",
            OptimizationResult(
                coordinates=ctx.molecule.coordinates,
                energy=-76.123456,
                converged=True,
            ),
        )
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert "-76.123456" in content

    def test_no_opt_annotation_when_absent(self, tmp_path):
        ctx = _make_context(tmp_path)
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert "Optimised energy" not in content

    def test_multipole_annotation_when_present(self, tmp_path):
        ctx = _make_context(tmp_path)
        ctx.set_artifact("esp_result", ESPGridResult(
            grid_coords=np.zeros((10, 3)), esp_values=np.zeros(10),
        ))
        ctx.set_artifact("multipole_frames", {
            "source": "dma", "dma_method": "MP2", "dma_basis_set": "6-311G**",
        })
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert "Multipole frames" in content
        assert "DMA" in content

    def test_torsion_annotation_when_present(self, tmp_path):
        ctx = _make_context(tmp_path)
        ctx.set_artifact("torsion_scans", {
            (0, 1): "scan_data",
            (2, 3): "scan_data",
        })
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert "Torsion scans: 2 bond(s)" in content

    def test_writes_fitted_parameter_sections(self, tmp_path):
        ctx = _make_context(tmp_path)
        ctx.set_artifact("multipoles", [
            MultipoleParam(
                atom_type=1,
                frame=[2],
                charge=-0.12345,
            )
        ])
        ctx.set_artifact("polarization_params", [
            PolarizeParam(atom_type=1, alpha=1.234, neighbors=[2])
        ])
        ctx.set_artifact("torsion_params", [
            TorsionParam(
                atom_types=(1, 2, 3, 4),
                folds=[TorsionFold(amplitude=0.5, phase=0.0, periodicity=3)],
            )
        ])
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert "multipole" in content
        assert "polarize" in content
        assert "torsion" in content

    def test_raises_on_no_molecule(self, tmp_path):
        config = PoltypeConfig()
        ctx = PipelineContext(config=config, work_dir=tmp_path)
        w = KeyFileWriter()
        with pytest.raises(ValueError, match="No molecule"):
            w.write(ctx)


# ===================================================================
# RunSummary tests
# ===================================================================


class TestRunSummary:
    """Tests for RunSummary dataclass."""

    def test_default_construction(self):
        s = RunSummary()
        assert s.molecule_name == "N/A"
        assert s.formula == "N/A"
        assert s.num_atoms == 0
        assert s.backend_name == "none"
        assert s.stage_results == {}
        assert s.output_files == []
        assert s.timestamp  # non-empty

    def test_from_context(self, tmp_path):
        ctx = _make_context(tmp_path)
        ctx.stage_results["input_preparation"] = StageResult(
            status=StageStatus.COMPLETED,
            message="OK",
            elapsed_seconds=0.5,
        )
        summary = RunSummary.from_context(ctx)
        assert summary.molecule_name == "ethanol"
        assert summary.num_atoms > 0
        assert summary.backend_name == "mock"
        assert "input_preparation" in summary.stage_results
        info = summary.stage_results["input_preparation"]
        assert info["status"] == "completed"
        assert info["elapsed_seconds"] == 0.5

    def test_from_context_with_output_files(self, tmp_path):
        ctx = _make_context(tmp_path)
        files = [Path("/tmp/a.xyz"), Path("/tmp/b.key")]
        summary = RunSummary.from_context(ctx, output_files=files)
        assert len(summary.output_files) == 2
        assert "/tmp/a.xyz" in summary.output_files

    def test_from_context_no_molecule(self, tmp_path):
        config = PoltypeConfig()
        ctx = PipelineContext(config=config, work_dir=tmp_path)
        summary = RunSummary.from_context(ctx)
        assert summary.molecule_name == "N/A"
        assert summary.num_atoms == 0

    def test_from_context_no_backend(self, tmp_path):
        config = PoltypeConfig()
        ctx = PipelineContext(
            config=config,
            molecule=_make_ethanol(),
            work_dir=tmp_path,
        )
        summary = RunSummary.from_context(ctx)
        assert summary.backend_name == "none"

    def test_timestamp_is_iso_format(self):
        s = RunSummary()
        # ISO format always contains 'T' separator
        assert "T" in s.timestamp


# ===================================================================
# SummaryWriter tests
# ===================================================================


class TestSummaryWriter:
    """Tests for the SummaryWriter."""

    def test_name(self):
        w = SummaryWriter()
        assert w.name == "summary_writer"

    def test_write_creates_file(self, tmp_path):
        ctx = _make_context(tmp_path)
        w = SummaryWriter()
        path = w.write(ctx)
        assert path.exists()
        assert path.name == "poltype_summary.txt"

    def test_write_custom_filename(self, tmp_path):
        ctx = _make_context(tmp_path)
        w = SummaryWriter(filename="report.txt")
        path = w.write(ctx)
        assert path.name == "report.txt"

    def test_content_has_molecule_info(self, tmp_path):
        ctx = _make_context(tmp_path)
        w = SummaryWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert "ethanol" in content
        assert "mock" in content

    def test_content_has_stage_results(self, tmp_path):
        ctx = _make_context(tmp_path)
        ctx.stage_results["input_preparation"] = StageResult(
            status=StageStatus.COMPLETED,
            message="OK",
            elapsed_seconds=1.23,
        )
        w = SummaryWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert "input_preparation" in content
        assert "COMPLETED" in content

    def test_content_has_output_files(self, tmp_path):
        ctx = _make_context(tmp_path)
        w = SummaryWriter()
        w.set_output_files([Path("/tmp/final.xyz"), Path("/tmp/final.key")])
        path = w.write(ctx)
        content = path.read_text()
        assert "final.xyz" in content
        assert "final.key" in content

    def test_set_output_files(self):
        w = SummaryWriter()
        w.set_output_files([Path("/a"), Path("/b")])
        assert len(w._extra_output_files) == 2


# ===================================================================
# FinalizationStage with writers tests
# ===================================================================


class TestFinalizationStageWithWriters:
    """Tests for the enhanced FinalizationStage."""

    def test_no_writers_backward_compatible(self, tmp_path):
        """Stage works without any writers (Phase 4 behaviour)."""
        ctx = _make_context(tmp_path)
        stage = FinalizationStage()
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        assert "final_summary" in result.artifacts

    def test_no_writers_no_output_files(self, tmp_path):
        """Without writers, output_files list is empty."""
        ctx = _make_context(tmp_path)
        stage = FinalizationStage()
        result = stage.execute(ctx)
        assert result.artifacts["output_files"] == []

    def test_with_writers_produces_files(self, tmp_path):
        """Writers produce actual files."""
        ctx = _make_context(tmp_path)
        writers = [XYZFileWriter(), KeyFileWriter()]
        stage = FinalizationStage(writers=writers)
        result = stage.execute(ctx)
        output_files = result.artifacts["output_files"]
        assert len(output_files) == 2
        assert all(p.exists() for p in output_files)

    def test_with_all_default_writers(self, tmp_path):
        """All three default writers produce files."""
        ctx = _make_context(tmp_path)
        writers = [XYZFileWriter(), KeyFileWriter(), SummaryWriter()]
        stage = FinalizationStage(writers=writers)
        result = stage.execute(ctx)
        output_files = result.artifacts["output_files"]
        assert len(output_files) == 3
        names = {p.name for p in output_files}
        assert "final.xyz" in names
        assert "final.key" in names
        assert "poltype_summary.txt" in names

    def test_writers_property(self):
        w = [XYZFileWriter()]
        stage = FinalizationStage(writers=w)
        assert len(stage.writers) == 1
        assert stage.writers[0].name == "xyz_writer"

    def test_writers_property_returns_copy(self):
        w = [XYZFileWriter()]
        stage = FinalizationStage(writers=w)
        stage.writers.append(KeyFileWriter())
        # Original should not change
        assert len(stage.writers) == 1

    def test_failing_writer_does_not_stop_others(self, tmp_path):
        """A writer that fails should not prevent others from running."""

        class FailingWriter(OutputWriter):
            @property
            def name(self):
                return "failing"

            def write(self, context):
                raise RuntimeError("Writer failed!")

        ctx = _make_context(tmp_path)
        writers = [FailingWriter(), XYZFileWriter()]
        stage = FinalizationStage(writers=writers)
        result = stage.execute(ctx)
        output_files = result.artifacts["output_files"]
        # Only XYZ should succeed
        assert len(output_files) == 1
        assert output_files[0].name == "final.xyz"

    def test_no_molecule_no_writer_output(self, tmp_path):
        """Without a molecule, writers are not invoked."""
        config = PoltypeConfig()
        ctx = PipelineContext(config=config, work_dir=tmp_path)
        writers = [XYZFileWriter(), KeyFileWriter()]
        stage = FinalizationStage(writers=writers)
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        assert result.artifacts["output_files"] == []

    def test_summary_includes_collected_artifacts(self, tmp_path):
        ctx = _make_context(tmp_path)
        ctx.set_artifact("opt_result", OptimizationResult(
            coordinates=ctx.molecule.coordinates, energy=-76.0, converged=True,
        ))
        stage = FinalizationStage()
        result = stage.execute(ctx)
        summary = result.artifacts["final_summary"]
        assert "opt_result" in summary["collected_artifacts"]
        assert summary["opt_energy"] == -76.0

    def test_result_message_includes_file_count(self, tmp_path):
        ctx = _make_context(tmp_path)
        writers = [XYZFileWriter()]
        stage = FinalizationStage(writers=writers)
        result = stage.execute(ctx)
        assert "1 file(s)" in result.message


# ===================================================================
# Dry-run mode tests
# ===================================================================


class TestDryRunMode:
    """Tests for dry-run mode across stages and pipeline."""

    def test_config_dry_run_default_false(self):
        config = PoltypeConfig()
        assert config.dry_run is False

    def test_config_dry_run_set_true(self):
        config = PoltypeConfig(dry_run=True)
        assert config.dry_run is True

    def test_geometry_opt_skips_in_dry_run(self, tmp_path):
        ctx = _make_context(tmp_path, dry_run=True)
        stage = GeometryOptimizationStage()
        assert stage.should_skip(ctx) is True

    def test_esp_skips_in_dry_run(self, tmp_path):
        ctx = _make_context(tmp_path, dry_run=True)
        stage = ESPFittingStage()
        assert stage.should_skip(ctx) is True

    def test_multipole_skips_in_dry_run(self, tmp_path):
        ctx = _make_context(tmp_path, dry_run=True)
        stage = MultipoleStage()
        assert stage.should_skip(ctx) is True

    def test_torsion_skips_in_dry_run(self, tmp_path):
        ctx = _make_context(tmp_path, dry_run=True)
        stage = TorsionFittingStage()
        assert stage.should_skip(ctx) is True

    def test_input_prep_does_not_skip_in_dry_run(self, tmp_path):
        """Input preparation should still run in dry-run mode."""
        ctx = _make_context(tmp_path, dry_run=True)
        stage = InputPreparationStage()
        assert stage.should_skip(ctx) is False

    def test_finalization_does_not_skip_in_dry_run(self, tmp_path):
        """Finalization should still run in dry-run mode."""
        ctx = _make_context(tmp_path, dry_run=True)
        stage = FinalizationStage()
        assert stage.should_skip(ctx) is False

    def test_geometry_opt_does_not_skip_normally(self, tmp_path):
        ctx = _make_context(tmp_path, dry_run=False)
        stage = GeometryOptimizationStage()
        assert stage.should_skip(ctx) is False

    def test_esp_does_not_skip_normally(self, tmp_path):
        ctx = _make_context(tmp_path, dry_run=False)
        stage = ESPFittingStage()
        assert stage.should_skip(ctx) is False

    def test_multipole_does_not_skip_normally(self, tmp_path):
        ctx = _make_context(tmp_path, dry_run=False)
        stage = MultipoleStage()
        assert stage.should_skip(ctx) is False

    def test_torsion_does_not_skip_normally(self, tmp_path):
        ctx = _make_context(tmp_path, dry_run=False)
        stage = TorsionFittingStage()
        assert stage.should_skip(ctx) is False

    def test_dry_run_full_pipeline(self, tmp_path):
        """Full pipeline in dry-run mode: only input_prep and finalization run."""
        mol = _make_ethanol()
        config = PoltypeConfig(dry_run=True)
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)

        # Input prep and finalization should have completed
        assert ctx.stage_results["input_preparation"].status is StageStatus.COMPLETED
        assert ctx.stage_results["finalization"].status is StageStatus.COMPLETED

        # QM stages and typing/database stages should have been skipped
        for name in ("geometry_optimization", "esp_fitting", "multipole",
                      "atom_typing", "database_match", "torsion_fitting"):
            assert ctx.stage_results[name].status is StageStatus.SKIPPED

    def test_dry_run_still_writes_output_files(self, tmp_path):
        """Dry-run should still produce output files from finalization."""
        mol = _make_ethanol()
        config = PoltypeConfig(dry_run=True)
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)

        output_files = ctx.artifacts.get("output_files", [])
        assert len(output_files) == 3  # xyz, key, summary
        assert all(p.exists() for p in output_files)


# ===================================================================
# Dry-run config loader integration tests
# ===================================================================


class TestDryRunConfigLoader:
    """Tests for dry-run in config loader."""

    def test_dry_run_from_ini(self, tmp_path):
        """dryrun keyword in poltype.ini sets dry_run=True."""
        from poltype.config.loader import load_config

        ini = tmp_path / "poltype.ini"
        ini.write_text("dryrun\n")
        cfg = load_config(ini)
        assert cfg.dry_run is True

    def test_dry_run_explicit_true(self, tmp_path):
        from poltype.config.loader import load_config

        ini = tmp_path / "poltype.ini"
        ini.write_text("dryrun = True\n")
        cfg = load_config(ini)
        assert cfg.dry_run is True

    def test_dry_run_explicit_false(self, tmp_path):
        from poltype.config.loader import load_config

        ini = tmp_path / "poltype.ini"
        ini.write_text("dryrun = False\n")
        cfg = load_config(ini)
        assert cfg.dry_run is False

    def test_dry_run_absent_defaults_false(self, tmp_path):
        from poltype.config.loader import load_config

        ini = tmp_path / "poltype.ini"
        ini.write_text("# no dryrun\n")
        cfg = load_config(ini)
        assert cfg.dry_run is False


# ===================================================================
# CLI --dry-run tests
# ===================================================================


class TestCLIDryRun:
    """Tests for the --dry-run CLI flag."""

    def test_parser_accepts_dry_run_flag(self):
        from poltype.__main__ import _build_parser

        parser = _build_parser()
        args = parser.parse_args(["--dry-run"])
        assert args.dry_run is True

    def test_parser_dry_run_default_false(self):
        from poltype.__main__ import _build_parser

        parser = _build_parser()
        args = parser.parse_args([])
        assert args.dry_run is False

    def test_parser_dry_run_with_other_flags(self):
        from poltype.__main__ import _build_parser

        parser = _build_parser()
        args = parser.parse_args(["--dry-run", "--verbose", "--no-checkpoint"])
        assert args.dry_run is True
        assert args.verbose is True
        assert args.no_checkpoint is True


# ===================================================================
# build_default_pipeline integration tests
# ===================================================================


class TestBuildDefaultPipeline:
    """Tests for the updated build_default_pipeline."""

    def test_returns_runner(self):
        runner = build_default_pipeline()
        assert isinstance(runner, PipelineRunner)

    def test_has_six_stages(self):
        runner = build_default_pipeline()
        assert len(runner._stages) == 10

    def test_finalization_has_writers(self):
        runner = build_default_pipeline()
        final_stage = runner._stages[-1]
        assert isinstance(final_stage, FinalizationStage)
        assert len(final_stage.writers) == 3

    def test_writer_names(self):
        runner = build_default_pipeline()
        final_stage = runner._stages[-1]
        writer_names = {w.name for w in final_stage.writers}
        assert writer_names == {"xyz_writer", "key_writer", "summary_writer"}

    def test_full_pipeline_produces_output_files(self, tmp_path):
        """Full pipeline (non-dry-run) with mock backend writes files."""
        mol = _make_ethanol()
        config = PoltypeConfig()
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)

        output_files = ctx.artifacts.get("output_files", [])
        assert len(output_files) == 3
        names = {p.name for p in output_files}
        assert "final.xyz" in names
        assert "final.key" in names
        assert "poltype_summary.txt" in names

    def test_full_pipeline_all_stages_complete(self, tmp_path):
        mol = _make_ethanol()
        config = PoltypeConfig()
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)

        for name in ("input_preparation", "geometry_optimization",
                      "esp_fitting", "multipole", "atom_typing",
                      "database_match", "torsion_fitting",
                      "finalization"):
            assert ctx.stage_results[name].status is StageStatus.COMPLETED


# ===================================================================
# Edge-case / regression tests
# ===================================================================


class TestEdgeCases:
    """Additional edge-case tests."""

    def test_xyz_writer_after_geometry_opt(self, tmp_path):
        """XYZ writer uses the latest molecule (post-optimization)."""
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        # Simulate geometry opt updating coordinates
        new_coords = mol.coordinates + 0.1
        new_mol = mol.with_updated_coordinates(new_coords)
        ctx.update_molecule(new_mol)

        w = XYZFileWriter()
        path = w.write(ctx)
        lines = path.read_text().strip().split("\n")
        # Check updated coords
        parts = lines[1].split()
        x = float(parts[2])
        assert abs(x - new_coords[0][0]) < 1e-5

    def test_key_writer_prm_reference(self, tmp_path):
        """Key file references the correct parameter file."""
        config = PoltypeConfig(prm_file_path="/custom/params.prm")
        ctx = PipelineContext(
            config=config,
            molecule=_make_ethanol(),
            work_dir=tmp_path,
        )
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert "/custom/params.prm" in content

    def test_summary_writer_with_no_stages(self, tmp_path):
        """Summary writer works even with no stage results."""
        ctx = _make_context(tmp_path)
        w = SummaryWriter()
        path = w.write(ctx)
        content = path.read_text()
        assert "Poltype2 Pipeline Run Summary" in content

    def test_dry_run_database_match_only_combined(self, tmp_path):
        """Both dry_run and database_match_only cause skipping."""
        config = PoltypeConfig(dry_run=True, database_match_only=True)
        ctx = PipelineContext(
            config=config,
            molecule=_make_ethanol(),
            backend=MockBackend(),
            work_dir=tmp_path,
        )
        stage = GeometryOptimizationStage()
        assert stage.should_skip(ctx) is True
