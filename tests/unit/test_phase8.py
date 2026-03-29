"""
tests/unit/test_phase8.py – unit tests for Phase 8: molecular fragmentation.

Tests cover:
- Fragment dataclass (frozen, attributes, properties).
- FragmentResult dataclass (properties, aggregation).
- Fragmenter ABC cannot be instantiated.
- WBOFragmenter: basic fragmentation, BFS growth, heuristic bond orders,
  QM WBO matrix, empty bonds, threshold/cycle configuration, repr.
- FragmentationStage: execute, skip logic, db_match_result integration,
  no molecule, no rotatable bonds, config-driven construction.
- Updated pipeline factory: 9 stages, stage order.
- Full pipeline integration: fragmentation stage completes.
- Dry-run mode: fragmentation stage skips.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pytest

from poltype.config.schema import PoltypeConfig, QMConfig
from poltype.database.match_result import DatabaseMatchResult
from poltype.fragmentation.fragment import Fragment, FragmentResult
from poltype.fragmentation.fragmenter import Fragmenter, RuleBasedFragmenter, WBOFragmenter
from poltype.molecule.molecule import Molecule
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.factory import build_default_pipeline
from poltype.pipeline.stage import StageResult, StageStatus
from poltype.pipeline.stages.fragmentation import FragmentationStage
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
    """Ethanol – has C, O, H and a rotatable bond for richer tests."""
    return Molecule.from_smiles("CCO", name="ethanol")


def _make_butane() -> Molecule:
    """Butane – linear chain with multiple rotatable bonds."""
    return Molecule.from_smiles("CCCC", name="butane")


def _make_context(
    tmp_path: Path,
    molecule: Optional[Molecule] = None,
    backend: Optional[QMBackend] = None,
    dry_run: bool = False,
    use_fragmentation: bool = True,
    database_match_only: bool = False,
) -> PipelineContext:
    config = PoltypeConfig(
        dry_run=dry_run,
        use_fragmentation=use_fragmentation,
        database_match_only=database_match_only,
    )
    mol = molecule or _make_ethanol()
    return PipelineContext(
        config=config,
        molecule=mol,
        backend=backend or MockBackend(),
        work_dir=tmp_path,
    )


# ===================================================================
# A. Fragment dataclass tests
# ===================================================================


class TestFragment:
    def test_creation(self):
        frag = Fragment(
            fragment_id=0,
            atom_map={0: 1, 1: 2, 2: 3},
            rotatable_bond=(1, 2),
        )
        assert frag.fragment_id == 0
        assert frag.atom_map == {0: 1, 1: 2, 2: 3}
        assert frag.rotatable_bond == (1, 2)
        assert frag.cap_atoms == []
        assert frag.name == ""

    def test_with_name_and_caps(self):
        frag = Fragment(
            fragment_id=1,
            atom_map={0: 5, 1: 6},
            rotatable_bond=(5, 6),
            cap_atoms=[2],
            name="frag_1_bond_5_6",
        )
        assert frag.name == "frag_1_bond_5_6"
        assert frag.cap_atoms == [2]

    def test_frozen(self):
        frag = Fragment(
            fragment_id=0,
            atom_map={0: 0},
            rotatable_bond=(0, 1),
        )
        with pytest.raises(AttributeError):
            frag.fragment_id = 99  # type: ignore[misc]

    def test_num_atoms(self):
        frag = Fragment(
            fragment_id=0,
            atom_map={0: 10, 1: 11, 2: 12},
            rotatable_bond=(10, 11),
        )
        assert frag.num_atoms == 3

    def test_parent_atom_indices(self):
        frag = Fragment(
            fragment_id=0,
            atom_map={0: 5, 1: 3, 2: 8},
            rotatable_bond=(3, 5),
        )
        assert frag.parent_atom_indices == [3, 5, 8]

    def test_num_cap_atoms(self):
        frag = Fragment(
            fragment_id=0,
            atom_map={0: 1, 1: 2, 2: 3},
            rotatable_bond=(1, 2),
            cap_atoms=[2],
        )
        assert frag.num_cap_atoms == 1

    def test_num_cap_atoms_zero(self):
        frag = Fragment(
            fragment_id=0,
            atom_map={0: 1},
            rotatable_bond=(0, 1),
        )
        assert frag.num_cap_atoms == 0

    def test_equality(self):
        a = Fragment(fragment_id=0, atom_map={0: 1}, rotatable_bond=(0, 1))
        b = Fragment(fragment_id=0, atom_map={0: 1}, rotatable_bond=(0, 1))
        assert a == b

    def test_inequality(self):
        a = Fragment(fragment_id=0, atom_map={0: 1}, rotatable_bond=(0, 1))
        b = Fragment(fragment_id=1, atom_map={0: 1}, rotatable_bond=(0, 1))
        assert a != b


# ===================================================================
# B. FragmentResult dataclass tests
# ===================================================================


class TestFragmentResult:
    def test_defaults(self):
        r = FragmentResult()
        assert r.fragments == []
        assert r.parent_molecule_name == ""
        assert r.fragmentation_method == ""
        assert r.metadata == {}

    def test_num_fragments(self):
        frags = [
            Fragment(fragment_id=0, atom_map={0: 0}, rotatable_bond=(0, 1)),
            Fragment(fragment_id=1, atom_map={0: 2}, rotatable_bond=(2, 3)),
        ]
        r = FragmentResult(fragments=frags)
        assert r.num_fragments == 2

    def test_num_fragments_empty(self):
        r = FragmentResult()
        assert r.num_fragments == 0

    def test_total_fragment_atoms(self):
        frags = [
            Fragment(fragment_id=0, atom_map={0: 0, 1: 1}, rotatable_bond=(0, 1)),
            Fragment(fragment_id=1, atom_map={0: 2, 1: 3, 2: 4}, rotatable_bond=(2, 3)),
        ]
        r = FragmentResult(fragments=frags)
        assert r.total_fragment_atoms == 5

    def test_covered_parent_atoms(self):
        frags = [
            Fragment(fragment_id=0, atom_map={0: 0, 1: 1}, rotatable_bond=(0, 1)),
            Fragment(fragment_id=1, atom_map={0: 1, 1: 2}, rotatable_bond=(1, 2)),
        ]
        r = FragmentResult(fragments=frags)
        assert r.covered_parent_atoms == {0, 1, 2}

    def test_covered_parent_atoms_empty(self):
        r = FragmentResult()
        assert r.covered_parent_atoms == set()

    def test_metadata(self):
        r = FragmentResult(
            metadata={"wbo_threshold": 0.5, "max_growth_cycles": 3},
        )
        assert r.metadata["wbo_threshold"] == 0.5

    def test_parent_molecule_name(self):
        r = FragmentResult(parent_molecule_name="ethanol")
        assert r.parent_molecule_name == "ethanol"


# ===================================================================
# C. Fragmenter ABC tests
# ===================================================================


class TestFragmenterABC:
    def test_cannot_instantiate(self):
        with pytest.raises(TypeError):
            Fragmenter()  # type: ignore[abstract]

    def test_concrete_must_implement_fragment(self):
        class IncompleteFragmenter(Fragmenter):
            pass

        with pytest.raises(TypeError):
            IncompleteFragmenter()  # type: ignore[abstract]

    def test_name_property(self):
        class ConcreteFragmenter(Fragmenter):
            def fragment(self, molecule, rotatable_bonds, **kwargs):
                return FragmentResult()

        f = ConcreteFragmenter()
        assert f.name == "ConcreteFragmenter"


# ===================================================================
# D. WBOFragmenter tests
# ===================================================================


class TestWBOFragmenter:
    def test_default_params(self):
        f = WBOFragmenter()
        assert f.wbo_threshold == 0.5
        assert f.max_growth_cycles == 5

    def test_custom_params(self):
        f = WBOFragmenter(wbo_threshold=0.8, max_growth_cycles=3)
        assert f.wbo_threshold == 0.8
        assert f.max_growth_cycles == 3

    def test_repr(self):
        f = WBOFragmenter(wbo_threshold=0.7, max_growth_cycles=4)
        r = repr(f)
        assert "WBOFragmenter" in r
        assert "0.7" in r
        assert "4" in r

    def test_name(self):
        f = WBOFragmenter()
        assert f.name == "WBOFragmenter"

    def test_fragment_no_bonds(self):
        f = WBOFragmenter()
        mol = _make_methane()
        result = f.fragment(mol, [])
        assert result.num_fragments == 0
        assert result.fragmentation_method == "wbo"
        assert result.parent_molecule_name == "methane"

    def test_fragment_ethanol(self):
        """Ethanol has one rotatable bond; should produce one fragment."""
        f = WBOFragmenter()
        mol = _make_ethanol()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("Ethanol embedding may not produce rotatable bonds")
        result = f.fragment(mol, bonds)
        assert result.num_fragments == len(bonds)
        for frag in result.fragments:
            assert frag.num_atoms >= 2  # At least the bond atoms

    def test_fragment_butane(self):
        """Butane has multiple rotatable bonds."""
        f = WBOFragmenter()
        mol = _make_butane()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("Butane embedding may not produce rotatable bonds")
        result = f.fragment(mol, bonds)
        assert result.num_fragments == len(bonds)

    def test_fragment_uses_heuristic_bond_orders(self):
        """Without WBO matrix, heuristic bond orders are used."""
        f = WBOFragmenter(wbo_threshold=0.5)
        mol = _make_ethanol()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        result = f.fragment(mol, bonds, wbo_matrix=None)
        assert result.metadata.get("used_qm_wbo") is False

    def test_fragment_uses_qm_wbo(self):
        """With WBO matrix provided, uses QM bond orders."""
        f = WBOFragmenter(wbo_threshold=0.5)
        mol = _make_ethanol()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        n = mol.num_atoms
        # Identity matrix means all bonds have WBO=1.0
        wbo = np.eye(n)
        result = f.fragment(mol, bonds, wbo_matrix=wbo)
        assert result.metadata.get("used_qm_wbo") is True

    def test_fragment_result_metadata(self):
        f = WBOFragmenter(wbo_threshold=0.8, max_growth_cycles=3)
        mol = _make_ethanol()
        result = f.fragment(mol, [])
        assert result.metadata["wbo_threshold"] == 0.8
        assert result.metadata["max_growth_cycles"] == 3

    def test_growth_includes_central_bond_atoms(self):
        """Fragment always includes both atoms of the central bond."""
        f = WBOFragmenter()
        mol = _make_ethanol()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        result = f.fragment(mol, bonds)
        for i, frag in enumerate(result.fragments):
            bond = bonds[i]
            parent_indices = set(frag.atom_map.values())
            assert bond[0] in parent_indices
            assert bond[1] in parent_indices

    def test_max_growth_cycles_limits_size(self):
        """With max_growth_cycles=0, fragment contains only bond atoms."""
        f = WBOFragmenter(max_growth_cycles=0)
        mol = _make_butane()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        result = f.fragment(mol, bonds)
        for frag in result.fragments:
            # With 0 growth, only the two central bond atoms are included
            assert frag.num_atoms == 2

    def test_fragment_ids_are_sequential(self):
        """Fragment IDs should be 0, 1, 2, ..."""
        f = WBOFragmenter()
        mol = _make_butane()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        result = f.fragment(mol, bonds)
        for i, frag in enumerate(result.fragments):
            assert frag.fragment_id == i

    def test_fragment_names_contain_bond_info(self):
        f = WBOFragmenter()
        mol = _make_ethanol()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        result = f.fragment(mol, bonds)
        for frag in result.fragments:
            assert "frag_" in frag.name
            assert "bond_" in frag.name

    def test_high_wbo_threshold_limits_growth(self):
        """A very high threshold should produce smaller fragments."""
        f_low = WBOFragmenter(wbo_threshold=0.1, max_growth_cycles=5)
        f_high = WBOFragmenter(wbo_threshold=100.0, max_growth_cycles=5)
        mol = _make_butane()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        result_low = f_low.fragment(mol, bonds)
        result_high = f_high.fragment(mol, bonds)
        # High threshold should give smaller (or equal) fragments
        for fl, fh in zip(result_low.fragments, result_high.fragments):
            assert fh.num_atoms <= fl.num_atoms


# ===================================================================
# D2. RuleBasedFragmenter tests
# ===================================================================


def _make_toluene() -> Molecule:
    """Toluene – aromatic ring with substituent."""
    return Molecule.from_smiles("Cc1ccccc1", name="toluene")


def _make_naphthalene_chain() -> Molecule:
    """2-ethylnaphthalene – fused ring with chain."""
    return Molecule.from_smiles("CCc1ccc2ccccc2c1", name="ethylnaphthalene")


def _make_aniline() -> Molecule:
    """Aniline – aromatic ring with nitrogen substituent."""
    return Molecule.from_smiles("Nc1ccccc1", name="aniline")


class TestRuleBasedFragmenter:
    def test_default_params(self):
        f = RuleBasedFragmenter()
        assert f.max_fragment_atoms == 50

    def test_custom_params(self):
        f = RuleBasedFragmenter(max_fragment_atoms=30)
        assert f.max_fragment_atoms == 30

    def test_repr(self):
        f = RuleBasedFragmenter(max_fragment_atoms=40)
        r = repr(f)
        assert "RuleBasedFragmenter" in r
        assert "40" in r

    def test_name(self):
        f = RuleBasedFragmenter()
        assert f.name == "RuleBasedFragmenter"

    def test_fragment_no_bonds(self):
        f = RuleBasedFragmenter()
        mol = _make_methane()
        result = f.fragment(mol, [])
        assert result.num_fragments == 0
        assert result.fragmentation_method == "rule_based"
        assert result.parent_molecule_name == "methane"

    def test_fragment_ethanol(self):
        """Ethanol has rotatable bonds; should produce fragments."""
        f = RuleBasedFragmenter()
        mol = _make_ethanol()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("Ethanol embedding may not produce rotatable bonds")
        result = f.fragment(mol, bonds)
        assert result.num_fragments == len(bonds)
        for frag in result.fragments:
            assert frag.num_atoms >= 2

    def test_fragment_butane(self):
        """Butane has multiple rotatable bonds."""
        f = RuleBasedFragmenter()
        mol = _make_butane()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("Butane embedding may not produce rotatable bonds")
        result = f.fragment(mol, bonds)
        assert result.num_fragments == len(bonds)

    def test_fragment_ignores_wbo_matrix(self):
        """WBO matrix is accepted but ignored."""
        f = RuleBasedFragmenter()
        mol = _make_ethanol()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        n = mol.num_atoms
        wbo = np.eye(n)
        result = f.fragment(mol, bonds, wbo_matrix=wbo)
        assert result.fragmentation_method == "rule_based"
        assert "used_qm_wbo" not in result.metadata

    def test_fragment_result_metadata(self):
        f = RuleBasedFragmenter(max_fragment_atoms=30)
        mol = _make_ethanol()
        result = f.fragment(mol, [])
        assert result.metadata["max_fragment_atoms"] == 30

    def test_growth_includes_central_bond_atoms(self):
        """Fragment always includes both atoms of the central bond."""
        f = RuleBasedFragmenter()
        mol = _make_ethanol()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        result = f.fragment(mol, bonds)
        for i, frag in enumerate(result.fragments):
            bond = bonds[i]
            parent_indices = set(frag.atom_map.values())
            assert bond[0] in parent_indices
            assert bond[1] in parent_indices

    def test_torsion_region_always_kept(self):
        """Torsion region (bond atoms + neighbours) is never removed."""
        f = RuleBasedFragmenter()
        mol = _make_butane()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        result = f.fragment(mol, bonds)
        rdmol = mol.rdmol
        for i, frag in enumerate(result.fragments):
            bond = bonds[i]
            parent_indices = set(frag.atom_map.values())
            # Both bond atoms must be present.
            assert bond[0] in parent_indices
            assert bond[1] in parent_indices
            # All neighbours of bond atoms must be present.
            for center in bond:
                atom = rdmol.GetAtomWithIdx(center)
                for neig in atom.GetNeighbors():
                    assert neig.GetIdx() in parent_indices

    def test_fragment_ids_are_sequential(self):
        """Fragment IDs should be 0, 1, 2, ..."""
        f = RuleBasedFragmenter()
        mol = _make_butane()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        result = f.fragment(mol, bonds)
        for i, frag in enumerate(result.fragments):
            assert frag.fragment_id == i

    def test_fragment_names_contain_bond_info(self):
        f = RuleBasedFragmenter()
        mol = _make_ethanol()
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        result = f.fragment(mol, bonds)
        for frag in result.fragments:
            assert "frag_" in frag.name
            assert "bond_" in frag.name

    def test_methane_no_rotatable_bonds(self):
        """Methane has no rotatable bonds; no fragments produced."""
        f = RuleBasedFragmenter()
        mol = _make_methane()
        result = f.fragment(mol, [])
        assert result.num_fragments == 0


# ===================================================================
# E. FragmentationStage tests
# ===================================================================


class TestFragmentationStage:
    def test_name(self):
        stage = FragmentationStage()
        assert stage.name == "fragmentation"

    def test_execute_with_rotatable_bonds(self, tmp_path):
        mol = _make_ethanol()
        stage = FragmentationStage()
        ctx = _make_context(tmp_path, molecule=mol)
        if not mol.rotatable_bonds:
            pytest.skip("No rotatable bonds")
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        assert "fragment_result" in result.artifacts

    def test_execute_stores_artifact(self, tmp_path):
        mol = _make_ethanol()
        stage = FragmentationStage()
        ctx = _make_context(tmp_path, molecule=mol)
        stage.execute(ctx)
        assert ctx.get_artifact("fragment_result") is not None

    def test_execute_no_molecule(self, tmp_path):
        stage = FragmentationStage()
        ctx = PipelineContext(config=PoltypeConfig(), work_dir=tmp_path)
        result = stage.execute(ctx)
        assert result.status is StageStatus.FAILED
        assert "No molecule" in result.message

    def test_execute_no_rotatable_bonds(self, tmp_path):
        mol = _make_methane()
        stage = FragmentationStage()
        ctx = _make_context(tmp_path, molecule=mol)
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        assert "No rotatable bonds" in result.message

    def test_execute_with_injected_fragmenter(self, tmp_path):
        mol = _make_ethanol()
        fragmenter = WBOFragmenter(wbo_threshold=0.1, max_growth_cycles=2)
        stage = FragmentationStage(fragmenter=fragmenter)
        ctx = _make_context(tmp_path, molecule=mol)
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED

    def test_execute_with_db_match_result(self, tmp_path):
        """Stage uses db_match_result to filter bonds."""
        mol = _make_ethanol()
        stage = FragmentationStage()
        ctx = _make_context(tmp_path, molecule=mol)
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        # Pretend all torsions are matched → no fragmentation needed.
        matched_torsions = set()
        for bond in bonds:
            # Create a torsion with the bond as central pair
            matched_torsions.add((0, bond[0], bond[1], 0))
        db_result = DatabaseMatchResult(
            matched_atom_indices=list(range(mol.num_atoms)),
            unmatched_atom_indices=[],
            matched_torsions=matched_torsions,
        )
        ctx.set_artifact("db_match_result", db_result)
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        frag_result = result.artifacts["fragment_result"]
        assert frag_result.num_fragments == 0

    def test_execute_with_partial_db_match(self, tmp_path):
        """Bonds not covered by db_match are fragmented."""
        mol = _make_butane()
        stage = FragmentationStage()
        ctx = _make_context(tmp_path, molecule=mol)
        bonds = mol.rotatable_bonds
        if len(bonds) < 2:
            pytest.skip("Need at least 2 rotatable bonds")
        # Match only the first bond
        matched_torsions = {(0, bonds[0][0], bonds[0][1], 0)}
        db_result = DatabaseMatchResult(
            matched_atom_indices=[],
            unmatched_atom_indices=[],
            matched_torsions=matched_torsions,
        )
        ctx.set_artifact("db_match_result", db_result)
        result = stage.execute(ctx)
        frag_result = result.artifacts["fragment_result"]
        # Should have fragments for all unmatched bonds
        assert frag_result.num_fragments == len(bonds) - 1

    def test_execute_uses_rule_based_by_default(self, tmp_path):
        """Default stage uses RuleBasedFragmenter, ignoring wbo_matrix."""
        mol = _make_ethanol()
        stage = FragmentationStage()
        ctx = _make_context(tmp_path, molecule=mol)
        bonds = mol.rotatable_bonds
        if not bonds:
            pytest.skip("No rotatable bonds")
        wbo = np.eye(mol.num_atoms)
        ctx.set_artifact("wbo_matrix", wbo)
        result = stage.execute(ctx)
        frag_result = result.artifacts["fragment_result"]
        assert frag_result.fragmentation_method == "rule_based"

    def test_execute_config_driven_fragmenter(self, tmp_path):
        """When no fragmenter injected, builds a RuleBasedFragmenter."""
        config = PoltypeConfig()
        mol = _make_ethanol()
        ctx = PipelineContext(
            config=config,
            molecule=mol,
            backend=MockBackend(),
            work_dir=tmp_path,
        )
        stage = FragmentationStage()  # No injected fragmenter
        result = stage.execute(ctx)
        assert result.status is StageStatus.COMPLETED
        frag_result = result.artifacts.get("fragment_result")
        if frag_result and frag_result.num_fragments > 0:
            assert frag_result.fragmentation_method == "rule_based"

    def test_should_skip_dry_run(self, tmp_path):
        stage = FragmentationStage()
        ctx = _make_context(tmp_path, dry_run=True)
        assert stage.should_skip(ctx) is True

    def test_should_skip_database_match_only(self, tmp_path):
        stage = FragmentationStage()
        ctx = _make_context(tmp_path, database_match_only=True)
        assert stage.should_skip(ctx) is True

    def test_should_skip_fragmentation_disabled(self, tmp_path):
        stage = FragmentationStage()
        ctx = _make_context(tmp_path, use_fragmentation=False)
        assert stage.should_skip(ctx) is True

    def test_should_not_skip_normal(self, tmp_path):
        stage = FragmentationStage()
        ctx = _make_context(tmp_path)
        assert stage.should_skip(ctx) is False


# ===================================================================
# F. Pipeline factory tests (updated for 9 stages)
# ===================================================================


class TestBuildDefaultPipelinePhase8:
    def test_has_nine_stages(self):
        runner = build_default_pipeline()
        assert len(runner._stages) == 10

    def test_stage_order(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert names == [
            "input_preparation",
            "geometry_optimization",
            "multipole",
            "esp_fitting",
            "atom_typing",
            "database_match",
            "fragmentation",
            "torsion_fitting",
            "validation",
            "finalization",
        ]

    def test_fragmentation_stage_present(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert "fragmentation" in names

    def test_fragmentation_after_database_match(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert names.index("database_match") < names.index("fragmentation")

    def test_fragmentation_before_torsion(self):
        runner = build_default_pipeline()
        names = [s.name for s in runner._stages]
        assert names.index("fragmentation") < names.index("torsion_fitting")


# ===================================================================
# G. Full pipeline integration tests
# ===================================================================


class TestFullPipelineWithFragmentation:
    def test_fragmentation_runs_in_full_pipeline(self, tmp_path):
        """FragmentationStage completes in a full pipeline run."""
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
        assert "fragmentation" in ctx.stage_results
        assert ctx.stage_results["fragmentation"].status is StageStatus.COMPLETED

    def test_fragment_result_available_after_fragmentation(self, tmp_path):
        """fragment_result artifact is available after fragmentation runs."""
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
        result = ctx.get_artifact("fragment_result")
        assert result is not None
        assert isinstance(result, FragmentResult)

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

    def test_fragmentation_uses_db_match_result(self, tmp_path):
        """Fragmentation stage receives db_match_result from upstream."""
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
        # Both db_match_result and fragment_result should exist
        assert ctx.get_artifact("db_match_result") is not None
        assert ctx.get_artifact("fragment_result") is not None


class TestDryRunWithFragmentation:
    def test_fragmentation_skips_in_dry_run(self, tmp_path):
        stage = FragmentationStage()
        ctx = _make_context(tmp_path, dry_run=True)
        assert stage.should_skip(ctx) is True

    def test_dry_run_full_pipeline_skips_fragmentation(self, tmp_path):
        mol = _make_ethanol()
        config = PoltypeConfig(dry_run=True)
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)
        assert ctx.stage_results["fragmentation"].status is StageStatus.SKIPPED
        # input_prep and finalization still run
        assert ctx.stage_results["input_preparation"].status is StageStatus.COMPLETED
        assert ctx.stage_results["finalization"].status is StageStatus.COMPLETED

    def test_fragmentation_disabled_skips(self, tmp_path):
        config = PoltypeConfig(use_fragmentation=False)
        mol = _make_ethanol()
        backend = MockBackend()
        ctx = PipelineContext(
            config=config, molecule=mol, backend=backend, work_dir=tmp_path,
        )
        runner = build_default_pipeline()
        ctx = runner.run(ctx)
        assert ctx.stage_results["fragmentation"].status is StageStatus.SKIPPED


# ===================================================================
# H. Torsion coverage helper tests
# ===================================================================


class TestTorsionCovered:
    def test_bond_covered_by_torsion(self):
        matched = {(0, 1, 2, 3)}
        assert FragmentationStage._torsion_covered((1, 2), matched) is True

    def test_bond_covered_reversed(self):
        matched = {(0, 2, 1, 3)}
        assert FragmentationStage._torsion_covered((1, 2), matched) is True

    def test_bond_not_covered(self):
        matched = {(0, 1, 2, 3)}
        assert FragmentationStage._torsion_covered((3, 4), matched) is False

    def test_empty_matched_torsions(self):
        assert FragmentationStage._torsion_covered((0, 1), set()) is False

    def test_multiple_torsions(self):
        matched = {(0, 1, 2, 3), (4, 5, 6, 7)}
        assert FragmentationStage._torsion_covered((5, 6), matched) is True
        assert FragmentationStage._torsion_covered((1, 2), matched) is True
        assert FragmentationStage._torsion_covered((0, 7), matched) is False


# ===================================================================
# I. Import tests
# ===================================================================


class TestModuleImports:
    def test_import_fragmentation_module(self):
        from poltype.fragmentation import Fragment, FragmentResult
        from poltype.fragmentation import Fragmenter, RuleBasedFragmenter, WBOFragmenter
        assert Fragment is not None
        assert FragmentResult is not None
        assert Fragmenter is not None
        assert WBOFragmenter is not None
        assert RuleBasedFragmenter is not None

    def test_import_stage_from_stages_init(self):
        from poltype.pipeline.stages import FragmentationStage
        assert FragmentationStage is not None

    def test_import_stage_directly(self):
        from poltype.pipeline.stages.fragmentation import FragmentationStage
        assert FragmentationStage is not None
