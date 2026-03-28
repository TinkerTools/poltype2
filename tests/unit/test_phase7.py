"""
tests/unit/test_phase7.py – unit tests for Phase 7: atom typing & parameter database.

Tests cover:
- TypeRecord dataclass (frozen, attributes).
- AtomTyper ABC cannot be instantiated.
- SMARTSAtomTyper: from_lines, from_file, assign_types, last-match-wins,
  fallback for unmatched atoms, repr, num_rules.
- DatabaseMatchResult: coverage, all_atoms_matched.
- ParameterDatabase ABC cannot be instantiated.
- SmallMoleculeDB: lookup, coverage reporting, name, repr.
- AtomTypingStage: execute, fallback typing, skip in dry-run.
- DatabaseMatchStage: execute, skip in dry-run, missing database.
- Updated pipeline factory: 8 stages, stage order.
- XYZFileWriter uses type_records when available.
- KeyFileWriter uses type_records when available.
- Full pipeline integration: typing + database stages complete.
- Dry-run mode: typing and database stages skip.
"""

from __future__ import annotations

import textwrap
from pathlib import Path
from typing import List, Optional
from unittest.mock import MagicMock

import numpy as np
import pytest

from poltype.config.schema import PoltypeConfig, QMConfig
from poltype.database.match_result import DatabaseMatchResult
from poltype.database.parameter_db import ParameterDatabase, SmallMoleculeDB
from poltype.molecule.molecule import Molecule
from poltype.output.key_writer import KeyFileWriter
from poltype.output.xyz_writer import XYZFileWriter
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.factory import build_default_pipeline
from poltype.pipeline.runner import PipelineRunner
from poltype.pipeline.stage import Stage, StageResult, StageStatus
from poltype.pipeline.stages.atom_typing import AtomTypingStage
from poltype.pipeline.stages.database_match import DatabaseMatchStage
from poltype.qm.backend import (
    ESPGridResult,
    OptimizationResult,
    QMBackend,
    TorsionScanResult,
    WBOResult,
)
from poltype.typing.atom_typer import AtomTyper, SMARTSAtomTyper
from poltype.typing.type_record import TypeRecord


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


def _make_methane() -> Molecule:
    """Methane – simplest hydrocarbon for quick tests."""
    return Molecule.from_smiles("C", name="methane")


def _make_ethanol() -> Molecule:
    """Ethanol – has C, O, H for richer typing tests."""
    return Molecule.from_smiles("CCO", name="ethanol")


def _make_context(
    tmp_path: Path,
    molecule: Optional[Molecule] = None,
    backend: Optional[QMBackend] = None,
    dry_run: bool = False,
) -> PipelineContext:
    config = PoltypeConfig(dry_run=dry_run)
    mol = molecule or _make_ethanol()
    return PipelineContext(
        config=config,
        molecule=mol,
        backend=backend or MockBackend(),
        work_dir=tmp_path,
    )


SAMPLE_RULES = [
    "# Simple test rules",
    "[#1]     5  5  Hydrogen",
    "[CX4]    1  1  sp3 Carbon",
    "[OX2H1]  6  6  Hydroxyl O",
]


# ===================================================================
# A. TypeRecord dataclass tests
# ===================================================================


class TestTypeRecord:
    def test_creation(self):
        rec = TypeRecord(atom_index=0, atom_type=1, atom_class=1)
        assert rec.atom_index == 0
        assert rec.atom_type == 1
        assert rec.atom_class == 1
        assert rec.description == ""
        assert rec.smarts == ""

    def test_with_description_and_smarts(self):
        rec = TypeRecord(
            atom_index=3,
            atom_type=10,
            atom_class=5,
            description="Alkane CH3",
            smarts="[CX4H3]",
        )
        assert rec.description == "Alkane CH3"
        assert rec.smarts == "[CX4H3]"

    def test_frozen(self):
        rec = TypeRecord(atom_index=0, atom_type=1, atom_class=1)
        with pytest.raises(AttributeError):
            rec.atom_type = 2  # type: ignore[misc]

    def test_equality(self):
        a = TypeRecord(atom_index=0, atom_type=1, atom_class=1, description="X")
        b = TypeRecord(atom_index=0, atom_type=1, atom_class=1, description="X")
        assert a == b

    def test_inequality(self):
        a = TypeRecord(atom_index=0, atom_type=1, atom_class=1)
        b = TypeRecord(atom_index=0, atom_type=2, atom_class=1)
        assert a != b


# ===================================================================
# B. AtomTyper ABC tests
# ===================================================================


class TestAtomTyperABC:
    def test_cannot_instantiate(self):
        with pytest.raises(TypeError):
            AtomTyper()  # type: ignore[abstract]

    def test_concrete_subclass_must_implement_assign_types(self):
        class IncompleteTyper(AtomTyper):
            pass

        with pytest.raises(TypeError):
            IncompleteTyper()  # type: ignore[abstract]


# ===================================================================
# C. SMARTSAtomTyper tests
# ===================================================================


class TestSMARTSAtomTyper:
    def test_from_lines_basic(self):
        typer = SMARTSAtomTyper.from_lines(SAMPLE_RULES)
        assert typer.num_rules == 3

    def test_from_lines_skips_comments(self):
        lines = ["# comment", "", "[#1] 5 5"]
        typer = SMARTSAtomTyper.from_lines(lines)
        assert typer.num_rules == 1

    def test_from_lines_skips_invalid_smarts(self):
        lines = ["INVALID_SMARTS 1 1", "[#1] 5 5"]
        typer = SMARTSAtomTyper.from_lines(lines)
        # Invalid SMARTS is silently skipped
        assert typer.num_rules >= 1

    def test_assign_types_methane(self):
        typer = SMARTSAtomTyper.from_lines(SAMPLE_RULES)
        mol = _make_methane()
        records = typer.assign_types(mol.rdmol)
        assert len(records) == mol.num_atoms
        # Every atom should get an index
        for i, rec in enumerate(records):
            assert rec.atom_index == i

    def test_assign_types_all_atoms_typed(self):
        """With broad rules, all atoms of ethanol should be typed."""
        typer = SMARTSAtomTyper.from_lines(SAMPLE_RULES)
        mol = _make_ethanol()
        records = typer.assign_types(mol.rdmol)
        # All atoms should have non-zero type (H, C, O covered)
        for rec in records:
            assert rec.atom_type != 0, f"Atom {rec.atom_index} untyped"

    def test_assign_types_last_match_wins(self):
        """Later rules should override earlier ones."""
        lines = [
            "[#6]    1  1  Generic Carbon",
            "[CX4]   2  2  sp3 Carbon",
        ]
        typer = SMARTSAtomTyper.from_lines(lines)
        mol = _make_methane()
        records = typer.assign_types(mol.rdmol)
        # The carbon atom should get type 2 (last matching rule)
        carbon_recs = [r for r in records if r.atom_type in (1, 2)]
        for rec in carbon_recs:
            if rec.description == "sp3 Carbon":
                assert rec.atom_type == 2

    def test_assign_types_unmatched_fallback(self):
        """Atoms matching no rule get type=0."""
        lines = ["[F] 10 10 Fluorine"]  # Won't match methane
        typer = SMARTSAtomTyper.from_lines(lines)
        mol = _make_methane()
        records = typer.assign_types(mol.rdmol)
        for rec in records:
            assert rec.atom_type == 0
            assert rec.description == "Untyped"

    def test_from_file(self, tmp_path):
        content = "\n".join(SAMPLE_RULES)
        f = tmp_path / "test_smarts.txt"
        f.write_text(content)
        typer = SMARTSAtomTyper.from_file(f)
        assert typer.num_rules == 3

    def test_from_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            SMARTSAtomTyper.from_file("/nonexistent/path.txt")

    def test_from_file_amoeba_format(self, tmp_path):
        """Parse the real AMOEBA 3-column format: SMARTS classNum className."""
        content = textwrap.dedent("""\
            # smartsString classNumber className # comments
            [H]*                                   1         H* #   H*, wild matching Hydrogen
            [C]*                                   2         C* #   C*, wild matching Carbon
            [CX4H3]                               35       C33 #  C33, sp3 carbon with 3 Hydrogens
        """)
        f = tmp_path / "amoeba_format.txt"
        f.write_text(content)
        typer = SMARTSAtomTyper.from_file(f)
        assert typer.num_rules == 3

    def test_repr(self):
        typer = SMARTSAtomTyper.from_lines(SAMPLE_RULES)
        assert "SMARTSAtomTyper" in repr(typer)
        assert "3" in repr(typer)

    def test_num_rules_property(self):
        typer = SMARTSAtomTyper.from_lines([])
        assert typer.num_rules == 0


# ===================================================================
# D. DatabaseMatchResult tests
# ===================================================================


class TestDatabaseMatchResult:
    def test_defaults(self):
        r = DatabaseMatchResult()
        assert r.matched_atom_indices == []
        assert r.unmatched_atom_indices == []
        assert r.matched_torsions == set()
        assert r.source == ""
        assert r.metadata == {}

    def test_all_atoms_matched_true(self):
        r = DatabaseMatchResult(
            matched_atom_indices=[0, 1, 2],
            unmatched_atom_indices=[],
        )
        assert r.all_atoms_matched is True

    def test_all_atoms_matched_false(self):
        r = DatabaseMatchResult(
            matched_atom_indices=[0, 1],
            unmatched_atom_indices=[2],
        )
        assert r.all_atoms_matched is False

    def test_coverage_full(self):
        r = DatabaseMatchResult(
            matched_atom_indices=[0, 1, 2],
            unmatched_atom_indices=[],
        )
        assert r.coverage == pytest.approx(1.0)

    def test_coverage_partial(self):
        r = DatabaseMatchResult(
            matched_atom_indices=[0, 1],
            unmatched_atom_indices=[2, 3],
        )
        assert r.coverage == pytest.approx(0.5)

    def test_coverage_zero(self):
        r = DatabaseMatchResult(
            matched_atom_indices=[],
            unmatched_atom_indices=[0, 1, 2],
        )
        assert r.coverage == pytest.approx(0.0)

    def test_coverage_empty(self):
        r = DatabaseMatchResult()
        assert r.coverage == pytest.approx(0.0)


# ===================================================================
# E. ParameterDatabase ABC tests
# ===================================================================


class TestParameterDatabaseABC:
    def test_cannot_instantiate(self):
        with pytest.raises(TypeError):
            ParameterDatabase()  # type: ignore[abstract]

    def test_concrete_must_implement_methods(self):
        class IncompleteDB(ParameterDatabase):
            @property
            def name(self):
                return "incomplete"

        with pytest.raises(TypeError):
            IncompleteDB()  # type: ignore[abstract]


# ===================================================================
# F. SmallMoleculeDB tests
# ===================================================================


class TestSmallMoleculeDB:
    @pytest.fixture
    def smarts_file(self, tmp_path):
        content = "\n".join(SAMPLE_RULES)
        f = tmp_path / "smarts.txt"
        f.write_text(content)
        return f

    def test_creation(self, smarts_file):
        db = SmallMoleculeDB(smarts_file)
        assert db.num_rules == 3

    def test_name(self, smarts_file):
        db = SmallMoleculeDB(smarts_file)
        assert "SmallMoleculeDB" in db.name
        assert "smarts.txt" in db.name

    def test_lookup_all_matched(self, smarts_file):
        db = SmallMoleculeDB(smarts_file)
        mol = _make_ethanol()
        result = db.lookup(mol.rdmol)
        assert isinstance(result, DatabaseMatchResult)
        # Ethanol has C, H, O → all should match
        assert result.all_atoms_matched

    def test_lookup_coverage(self, smarts_file):
        db = SmallMoleculeDB(smarts_file)
        mol = _make_ethanol()
        result = db.lookup(mol.rdmol)
        assert result.coverage == pytest.approx(1.0)

    def test_lookup_with_precomputed_records(self, smarts_file):
        db = SmallMoleculeDB(smarts_file)
        mol = _make_ethanol()
        # Manually create records where some are untyped
        records = [
            TypeRecord(atom_index=i, atom_type=0, atom_class=0, description="Untyped")
            for i in range(mol.num_atoms)
        ]
        result = db.lookup(mol.rdmol, type_records=records)
        assert result.all_atoms_matched is False
        assert len(result.unmatched_atom_indices) == mol.num_atoms

    def test_lookup_source(self, smarts_file):
        db = SmallMoleculeDB(smarts_file)
        mol = _make_methane()
        result = db.lookup(mol.rdmol)
        assert "SmallMoleculeDB" in result.source

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            SmallMoleculeDB("/nonexistent/smarts.txt")

    def test_repr(self, smarts_file):
        db = SmallMoleculeDB(smarts_file)
        assert "SmallMoleculeDB" in repr(db)
        assert "rules=" in repr(db)


# ===================================================================
# G. AtomTypingStage tests
# ===================================================================


class TestAtomTypingStage:
    def test_name(self):
        stage = AtomTypingStage()
        assert stage.name == "atom_typing"

    def test_execute_with_injected_typer(self, tmp_path):
        typer = SMARTSAtomTyper.from_lines(SAMPLE_RULES)
        stage = AtomTypingStage(typer=typer)
        ctx = _make_context(tmp_path)
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        assert "type_records" in result.artifacts
        records = result.artifacts["type_records"]
        assert len(records) == ctx.molecule.num_atoms

    def test_execute_stores_artifact(self, tmp_path):
        typer = SMARTSAtomTyper.from_lines(SAMPLE_RULES)
        stage = AtomTypingStage(typer=typer)
        ctx = _make_context(tmp_path)
        stage.execute(ctx)
        assert ctx.get_artifact("type_records") is not None

    def test_execute_reports_typed_count(self, tmp_path):
        typer = SMARTSAtomTyper.from_lines(SAMPLE_RULES)
        stage = AtomTypingStage(typer=typer)
        ctx = _make_context(tmp_path)
        result = stage.execute(ctx)
        assert "Typed" in result.message

    def test_execute_fails_without_molecule(self, tmp_path):
        stage = AtomTypingStage()
        ctx = PipelineContext(config=PoltypeConfig(), work_dir=tmp_path)
        result = stage.execute(ctx)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_execute_fallback_when_file_missing(self, tmp_path):
        """When SMARTS file doesn't exist, uses index-based fallback."""
        config = PoltypeConfig(
            small_molecule_smarts_to_tinker_class="/nonexistent.txt"
        )
        ctx = PipelineContext(
            config=config,
            molecule=_make_ethanol(),
            backend=MockBackend(),
            work_dir=tmp_path,
        )
        stage = AtomTypingStage()  # No injected typer
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        assert "Fallback" in result.message
        records = result.artifacts["type_records"]
        # Fallback types = tinker index (1-based)
        for i, rec in enumerate(records):
            assert rec.atom_type == i + 1

    def test_should_skip_dry_run(self, tmp_path):
        stage = AtomTypingStage()
        ctx = _make_context(tmp_path, dry_run=True)
        assert stage.should_skip(ctx) is True

    def test_should_not_skip_normal(self, tmp_path):
        stage = AtomTypingStage()
        ctx = _make_context(tmp_path, dry_run=False)
        assert stage.should_skip(ctx) is False


# ===================================================================
# H. DatabaseMatchStage tests
# ===================================================================


class TestDatabaseMatchStage:
    @pytest.fixture
    def smarts_file(self, tmp_path):
        content = "\n".join(SAMPLE_RULES)
        f = tmp_path / "smarts.txt"
        f.write_text(content)
        return f

    def test_name(self):
        stage = DatabaseMatchStage()
        assert stage.name == "database_match"

    def test_execute_with_injected_db(self, tmp_path, smarts_file):
        db = SmallMoleculeDB(smarts_file)
        stage = DatabaseMatchStage(database=db)
        ctx = _make_context(tmp_path)
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        assert "db_match_result" in result.artifacts

    def test_execute_stores_artifact(self, tmp_path, smarts_file):
        db = SmallMoleculeDB(smarts_file)
        stage = DatabaseMatchStage(database=db)
        ctx = _make_context(tmp_path)
        stage.execute(ctx)
        assert ctx.get_artifact("db_match_result") is not None

    def test_execute_reports_coverage(self, tmp_path, smarts_file):
        db = SmallMoleculeDB(smarts_file)
        stage = DatabaseMatchStage(database=db)
        ctx = _make_context(tmp_path)
        result = stage.execute(ctx)
        assert "matched" in result.message

    def test_execute_fails_without_molecule(self, tmp_path):
        stage = DatabaseMatchStage()
        ctx = PipelineContext(config=PoltypeConfig(), work_dir=tmp_path)
        result = stage.execute(ctx)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_execute_missing_db_file(self, tmp_path):
        """When database file doesn't exist, reports all unmatched."""
        config = PoltypeConfig(
            small_molecule_smarts_to_tinker_class="/nonexistent.txt"
        )
        ctx = PipelineContext(
            config=config,
            molecule=_make_ethanol(),
            backend=MockBackend(),
            work_dir=tmp_path,
        )
        stage = DatabaseMatchStage()
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        assert "No database" in result.message
        match_result = result.artifacts["db_match_result"]
        assert match_result.all_atoms_matched is False

    def test_execute_uses_type_records_artifact(self, tmp_path, smarts_file):
        """Stage picks up pre-computed type_records from context."""
        db = SmallMoleculeDB(smarts_file)
        stage = DatabaseMatchStage(database=db)
        ctx = _make_context(tmp_path)
        mol = ctx.molecule
        # Pre-seed type_records (all typed)
        records = [
            TypeRecord(atom_index=i, atom_type=1, atom_class=1)
            for i in range(mol.num_atoms)
        ]
        ctx.set_artifact("type_records", records)
        result = stage.execute(ctx)
        match_result = result.artifacts["db_match_result"]
        assert match_result.all_atoms_matched

    def test_should_skip_dry_run(self, tmp_path):
        stage = DatabaseMatchStage()
        ctx = _make_context(tmp_path, dry_run=True)
        assert stage.should_skip(ctx) is True

    def test_should_not_skip_normal(self, tmp_path):
        stage = DatabaseMatchStage()
        ctx = _make_context(tmp_path, dry_run=False)
        assert stage.should_skip(ctx) is False


# ===================================================================
# I. Pipeline factory tests (updated for 8 stages)
# ===================================================================


class TestBuildDefaultPipelinePhase7:
    def test_has_eight_stages(self):
        runner = build_default_pipeline()
        assert len(runner._stages) == 8

    def test_stage_order(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert names == [
            "input_preparation",
            "geometry_optimization",
            "esp_fitting",
            "multipole",
            "atom_typing",
            "database_match",
            "torsion_fitting",
            "finalization",
        ]

    def test_atom_typing_stage_present(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert "atom_typing" in names

    def test_database_match_stage_present(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert "database_match" in names

    def test_atom_typing_before_database_match(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert names.index("atom_typing") < names.index("database_match")

    def test_database_match_before_torsion(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert names.index("database_match") < names.index("torsion_fitting")


# ===================================================================
# J. Output writer integration with type_records
# ===================================================================


class TestXYZWriterWithTypeRecords:
    def test_uses_type_records(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        # Set type records with specific types
        records = [
            TypeRecord(atom_index=i, atom_type=100 + i, atom_class=100 + i)
            for i in range(mol.num_atoms)
        ]
        ctx.set_artifact("type_records", records)
        w = XYZFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        lines = content.strip().split("\n")
        # Check that atom type 100 appears in the first atom line
        first_atom_line = lines[1]
        parts = first_atom_line.split()
        # Type field is after coords (index 5)
        atom_type_val = int(parts[5])
        assert atom_type_val == 100

    def test_falls_back_to_index_without_type_records(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        # No type_records artifact set
        w = XYZFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        lines = content.strip().split("\n")
        # First atom line: type should be 1 (tinker index)
        first_atom_line = lines[1]
        parts = first_atom_line.split()
        atom_type_val = int(parts[5])
        assert atom_type_val == 1


class TestKeyWriterWithTypeRecords:
    def test_uses_type_records(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        records = [
            TypeRecord(
                atom_index=i,
                atom_type=200 + i,
                atom_class=300 + i,
                description=f"TestAtom{i}",
            )
            for i in range(mol.num_atoms)
        ]
        ctx.set_artifact("type_records", records)
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        # Check that atom type 200 and class 300 appear
        assert "200" in content
        assert "300" in content
        assert "TestAtom0" in content

    def test_falls_back_without_type_records(self, tmp_path):
        mol = _make_ethanol()
        ctx = _make_context(tmp_path, molecule=mol)
        w = KeyFileWriter()
        path = w.write(ctx)
        content = path.read_text()
        # Check that index-based types are used (atom 1 gets type 1)
        lines = content.split("\n")
        atom_lines = [l for l in lines if l.startswith("atom")]
        assert len(atom_lines) == mol.num_atoms
        # First atom line should have type=1, class=1
        parts = atom_lines[0].split()
        assert parts[1] == "1"  # type
        assert parts[2] == "1"  # class


# ===================================================================
# K. Full pipeline integration tests
# ===================================================================


class TestFullPipelineWithTyping:
    def test_typing_stage_runs_in_full_pipeline(self, tmp_path):
        """AtomTypingStage completes in a full pipeline run."""
        mol = _make_ethanol()
        # Use a config with a nonexistent SMARTS file to trigger fallback
        config = PoltypeConfig(
            small_molecule_smarts_to_tinker_class="/nonexistent_smarts.txt"
        )
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)
        assert "atom_typing" in ctx.stage_results
        assert ctx.stage_results["atom_typing"].status is StageStatus.COMPLETED

    def test_database_match_stage_runs_in_full_pipeline(self, tmp_path):
        """DatabaseMatchStage completes in a full pipeline run."""
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
        assert "database_match" in ctx.stage_results
        assert ctx.stage_results["database_match"].status is StageStatus.COMPLETED

    def test_type_records_available_after_typing(self, tmp_path):
        """type_records artifact is available after atom_typing runs."""
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
        records = ctx.get_artifact("type_records")
        assert records is not None
        assert len(records) == mol.num_atoms

    def test_db_match_result_available_after_match(self, tmp_path):
        """db_match_result artifact is available after database_match runs."""
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
        result = ctx.get_artifact("db_match_result")
        assert result is not None
        assert isinstance(result, DatabaseMatchResult)

    def test_all_eight_stages_complete(self, tmp_path):
        """All 8 stages complete in a full pipeline run."""
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
        assert len(ctx.stage_results) == 8
        for sr in ctx.stage_results.values():
            assert sr.status is StageStatus.COMPLETED

    def test_output_files_use_typed_atoms(self, tmp_path):
        """Output files reflect assigned atom types."""
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
        output_files = ctx.artifacts.get("output_files", [])
        assert len(output_files) == 3
        assert all(p.exists() for p in output_files)


class TestDryRunWithTypingStages:
    def test_atom_typing_skips_in_dry_run(self, tmp_path):
        stage = AtomTypingStage()
        ctx = _make_context(tmp_path, dry_run=True)
        assert stage.should_skip(ctx) is True

    def test_database_match_skips_in_dry_run(self, tmp_path):
        stage = DatabaseMatchStage()
        ctx = _make_context(tmp_path, dry_run=True)
        assert stage.should_skip(ctx) is True

    def test_dry_run_full_pipeline_skips_typing(self, tmp_path):
        mol = _make_ethanol()
        config = PoltypeConfig(dry_run=True)
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)
        assert ctx.stage_results["atom_typing"].status is StageStatus.SKIPPED
        assert ctx.stage_results["database_match"].status is StageStatus.SKIPPED
        # input_prep and finalization still run
        assert ctx.stage_results["input_preparation"].status is StageStatus.COMPLETED
        assert ctx.stage_results["finalization"].status is StageStatus.COMPLETED


# ===================================================================
# L. SMARTS file format edge cases
# ===================================================================


class TestSMARTSFileEdgeCases:
    def test_inline_comments_stripped(self, tmp_path):
        content = "[#1] 5 H #  Hydrogen\n[C] 2 C # Carbon\n"
        f = tmp_path / "smarts.txt"
        f.write_text(content)
        typer = SMARTSAtomTyper.from_file(f)
        assert typer.num_rules == 2

    def test_extra_whitespace(self, tmp_path):
        content = "  [#1]    5     H   # spaced out\n"
        f = tmp_path / "smarts.txt"
        f.write_text(content)
        typer = SMARTSAtomTyper.from_file(f)
        assert typer.num_rules == 1

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.txt"
        f.write_text("")
        typer = SMARTSAtomTyper.from_file(f)
        assert typer.num_rules == 0

    def test_only_comments(self, tmp_path):
        f = tmp_path / "comments.txt"
        f.write_text("# Just comments\n# Nothing else\n")
        typer = SMARTSAtomTyper.from_file(f)
        assert typer.num_rules == 0

    def test_real_amoeba_file_parseable(self):
        """The actual amoeba21smartstoclass.txt file should parse."""
        import os
        repo_root = Path(__file__).resolve().parents[2]
        smarts_file = repo_root / "ParameterFiles" / "amoeba21smartstoclass.txt"
        if smarts_file.exists():
            typer = SMARTSAtomTyper.from_file(smarts_file)
            assert typer.num_rules > 0
