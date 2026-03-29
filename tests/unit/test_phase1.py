"""
tests/unit/test_phase1.py – unit tests for Phase 1: Configuration system
and Molecule representation (data layer).

Tests cover:
- PoltypeConfig schema, defaults, and sub-config wiring.
- QMConfig and ResourceConfig defaults and field constraints.
- INI parser edge cases (bare keywords, comments, type coercion).
- load_config() integration with config file on disk.
- Molecule construction from SMILES, files, and RDKit mol objects.
- Molecule property accessors (num_atoms, formula, charge, coordinates).
- Index conversion between Tinker (1-based) and RDKit (0-based).
- Symmetry class computation and rotatable bond detection.
- Immutable coordinate update returning a new Molecule.
- Molecule equality and repr.
- Config ↔ Molecule integration (structure path → Molecule).
"""

from __future__ import annotations

import os
import textwrap
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from poltype.config.loader import _parse_ini_lines, load_config
from poltype.config.schema import PoltypeConfig, QMConfig, ResourceConfig
from poltype.molecule.molecule import Molecule


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def methane() -> Molecule:
    """CH4 – 5 atoms, no rotatable bonds."""
    return Molecule.from_smiles("C", name="methane", embed_3d=True)


@pytest.fixture
def ethanol() -> Molecule:
    """Ethanol with 3D coords."""
    return Molecule.from_smiles("CCO", name="ethanol", embed_3d=True)


@pytest.fixture
def default_config() -> PoltypeConfig:
    return PoltypeConfig()


@pytest.fixture
def sdf_file(tmp_path) -> Path:
    """Write a minimal methane SDF to a temp file."""
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    sdf = tmp_path / "methane.sdf"
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()
    return sdf


# ===================================================================
# A. PoltypeConfig schema tests
# ===================================================================


class TestPoltypeConfigSchema:
    """Validate the top-level config schema and its sub-configs."""

    def test_instantiation_with_no_args(self):
        cfg = PoltypeConfig()
        assert cfg is not None

    def test_qm_sub_config_is_qm_config(self):
        cfg = PoltypeConfig()
        assert isinstance(cfg.qm, QMConfig)

    def test_resource_sub_config_is_resource_config(self):
        cfg = PoltypeConfig()
        assert isinstance(cfg.resources, ResourceConfig)

    def test_force_field_default_amoeba(self):
        cfg = PoltypeConfig()
        assert cfg.force_field == "AMOEBA"

    def test_structure_default_none(self):
        cfg = PoltypeConfig()
        assert cfg.structure is None

    def test_prm_file_path_is_absolute(self):
        cfg = PoltypeConfig()
        assert os.path.isabs(cfg.prm_file_path)

    def test_database_match_only_default_false(self):
        cfg = PoltypeConfig()
        assert cfg.database_match_only is False

    def test_do_torsion_fit_default_true(self):
        cfg = PoltypeConfig()
        assert cfg.do_torsion_fit is True

    def test_total_charge_default_none(self):
        cfg = PoltypeConfig()
        assert cfg.total_charge is None


class TestQMConfigDefaults:
    """QM sub-config field defaults."""

    def test_opt_method(self):
        assert QMConfig().opt_method == "MP2"

    def test_opt_basis_set(self):
        assert QMConfig().opt_basis_set == "6-31G*"

    def test_backend(self):
        assert QMConfig().backend == "psi4"

    def test_use_gaus_false(self):
        assert QMConfig().use_gaus is False

    def test_dont_use_pyscf_false(self):
        assert QMConfig().dont_use_pyscf is False

    def test_opt_loose_true(self):
        assert QMConfig().opt_loose is True


class TestResourceConfigDefaults:
    """Resource sub-config field defaults."""

    def test_num_proc_none(self):
        assert ResourceConfig().num_proc is None

    def test_max_mem_gb(self):
        assert ResourceConfig().max_mem_gb == 16

    def test_consumption_ratio(self):
        assert ResourceConfig().consumption_ratio == 0.8


# ===================================================================
# B. INI parser tests
# ===================================================================


class TestParseIniLines:
    """Low-level _parse_ini_lines edge cases."""

    def test_simple_key_value(self):
        assert _parse_ini_lines(["k = v\n"])["k"] == "v"

    def test_comment_lines_ignored(self):
        result = _parse_ini_lines(["# a comment\n", "key = val\n"])
        assert len(result) == 1
        assert result["key"] == "val"

    def test_inline_comment_stripped(self):
        result = _parse_ini_lines(["key = val  # inline\n"])
        assert result["key"] == "val"

    def test_bare_keyword_maps_to_none(self):
        result = _parse_ini_lines(["debugmode\n"])
        assert "debugmode" in result
        assert result["debugmode"] is None

    def test_none_value_excluded(self):
        result = _parse_ini_lines(["key = None\n"])
        assert "key" not in result

    def test_empty_lines_ignored(self):
        result = _parse_ini_lines(["\n", "  \n", "x = 1\n"])
        assert len(result) == 1

    def test_keys_lowercased(self):
        result = _parse_ini_lines(["OptMethod = MP2\n"])
        assert "optmethod" in result

    def test_whitespace_around_value(self):
        result = _parse_ini_lines(["key =  spaced  \n"])
        assert result["key"] == "spaced"

    def test_multiple_assignments(self):
        lines = ["a = 1\n", "b = 2\n", "c = 3\n"]
        result = _parse_ini_lines(lines)
        assert len(result) == 3
        assert result["a"] == "1"
        assert result["b"] == "2"
        assert result["c"] == "3"


# ===================================================================
# C. load_config integration
# ===================================================================


class TestLoadConfigIntegration:
    """load_config from a real file on disk."""

    def _write_ini(self, tmp_path: Path, content: str) -> Path:
        ini = tmp_path / "poltype.ini"
        ini.write_text(textwrap.dedent(content))
        return ini

    def test_empty_file_returns_default_config(self, tmp_path):
        ini = self._write_ini(tmp_path, "")
        cfg = load_config(ini)
        assert isinstance(cfg, PoltypeConfig)

    def test_opt_method_override(self, tmp_path):
        ini = self._write_ini(tmp_path, "optmethod = wB97X-D\n")
        cfg = load_config(ini)
        assert cfg.qm.opt_method == "wB97X-D"

    def test_numproc_parsed_as_int(self, tmp_path):
        ini = self._write_ini(tmp_path, "numproc = 4\n")
        cfg = load_config(ini)
        assert cfg.resources.num_proc == 4
        assert isinstance(cfg.resources.num_proc, int)

    def test_totalcharge_parsed_as_int(self, tmp_path):
        ini = self._write_ini(tmp_path, "totalcharge = -1\n")
        cfg = load_config(ini)
        assert cfg.total_charge == -1

    def test_bool_bare_keyword(self, tmp_path):
        ini = self._write_ini(tmp_path, "debugmode\n")
        cfg = load_config(ini)
        assert cfg.debug_mode is True

    def test_use_gaus_flag(self, tmp_path):
        ini = self._write_ini(tmp_path, "use_gaus = True\n")
        cfg = load_config(ini)
        assert cfg.qm.use_gaus is True

    def test_list_field_comma_separated(self, tmp_path):
        ini = self._write_ini(tmp_path, "vdwprmtypestofit = S, T, R\n")
        cfg = load_config(ini)
        assert cfg.vdw_prm_types_to_fit == ["S", "T", "R"]

    def test_int_list_field(self, tmp_path):
        ini = self._write_ini(tmp_path, "onlyvdwatomlist = 1, 3, 5\n")
        cfg = load_config(ini)
        assert cfg.only_vdw_atom_list == [1, 3, 5]

    def test_unknown_keys_silently_ignored(self, tmp_path):
        ini = self._write_ini(
            tmp_path, "totally_unknown_key = something\n"
        )
        cfg = load_config(ini)
        assert isinstance(cfg, PoltypeConfig)

    def test_file_not_found_raises(self):
        with pytest.raises(FileNotFoundError):
            load_config("/nonexistent/poltype.ini")

    def test_mixed_comments_and_values(self, tmp_path):
        ini = self._write_ini(
            tmp_path,
            """\
            # QM settings
            optmethod = wB97X-D
            # resource settings
            numproc = 8
            """,
        )
        cfg = load_config(ini)
        assert cfg.qm.opt_method == "wB97X-D"
        assert cfg.resources.num_proc == 8

    def test_optloose_sets_convergence(self, tmp_path):
        ini = self._write_ini(tmp_path, "optloose = True\n")
        cfg = load_config(ini)
        assert cfg.qm.opt_loose is True
        assert cfg.qm.opt_convergence == "LOOSE"


# ===================================================================
# D. Molecule construction tests
# ===================================================================


class TestMoleculeConstruction:
    def test_from_smiles(self):
        mol = Molecule.from_smiles("C", embed_3d=False)
        assert mol is not None

    def test_from_smiles_with_name(self):
        mol = Molecule.from_smiles("C", name="methane", embed_3d=False)
        assert mol.name == "methane"

    def test_from_smiles_with_embed(self, methane):
        assert methane.has_conformer

    def test_from_rdmol(self):
        rdmol = Chem.MolFromSmiles("C")
        rdmol = Chem.AddHs(rdmol)
        mol = Molecule(rdmol, name="methane_direct")
        assert mol.num_atoms == 5

    def test_from_file_sdf(self, sdf_file):
        mol = Molecule.from_file(sdf_file)
        assert mol is not None
        assert mol.num_atoms > 0

    def test_from_file_default_name(self, sdf_file):
        mol = Molecule.from_file(sdf_file)
        assert mol.name == "methane"

    def test_from_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            Molecule.from_file("/nonexistent/molecule.sdf")

    def test_from_file_unsupported_format(self, tmp_path):
        p = tmp_path / "mol.xyz"
        p.write_text("dummy")
        with pytest.raises(ValueError, match="Unsupported file format"):
            Molecule.from_file(p)

    def test_none_rdmol_raises(self):
        with pytest.raises(ValueError):
            Molecule(None)


# ===================================================================
# E. Molecule properties
# ===================================================================


class TestMoleculeProperties:
    def test_num_atoms(self, methane):
        assert methane.num_atoms == 5  # C + 4H

    def test_num_heavy_atoms(self, methane):
        assert methane.num_heavy_atoms == 1

    def test_formal_charge_neutral(self, methane):
        assert methane.formal_charge == 0

    def test_formal_charge_cation(self):
        mol = Molecule.from_smiles("[NH4+]", embed_3d=False)
        assert mol.formal_charge == 1

    def test_formula(self, methane):
        assert methane.formula == "CH4"

    def test_has_conformer(self, methane):
        assert methane.has_conformer is True

    def test_no_conformer(self):
        mol = Molecule.from_smiles("C", embed_3d=False)
        assert mol.has_conformer is False

    def test_coordinates_shape(self, methane):
        assert methane.coordinates.shape == (5, 3)

    def test_coordinates_no_conformer_raises(self):
        mol = Molecule.from_smiles("C", embed_3d=False)
        with pytest.raises(RuntimeError, match="no conformer"):
            _ = mol.coordinates

    def test_rdmol_is_chem_mol(self, methane):
        assert isinstance(methane.rdmol, Chem.Mol)

    def test_work_dir_is_absolute(self, methane):
        assert methane.work_dir.is_absolute()


# ===================================================================
# F. Index conversion
# ===================================================================


class TestIndexConversion:
    def test_to_tinker_idx(self):
        assert Molecule.to_tinker_idx(0) == 1
        assert Molecule.to_tinker_idx(4) == 5

    def test_to_rdkit_idx(self):
        assert Molecule.to_rdkit_idx(1) == 0
        assert Molecule.to_rdkit_idx(5) == 4

    def test_roundtrip(self):
        for i in range(10):
            assert Molecule.to_rdkit_idx(Molecule.to_tinker_idx(i)) == i


# ===================================================================
# G. Symmetry classes
# ===================================================================


class TestSymmetryClasses:
    def test_returns_dict(self, methane):
        sym = methane.symmetry_classes
        assert isinstance(sym, dict)

    def test_all_atoms_present(self, methane):
        sym = methane.symmetry_classes
        assert set(sym.keys()) == set(range(methane.num_atoms))

    def test_equivalent_hydrogens(self, methane):
        sym = methane.symmetry_classes
        h_ranks = {sym[i] for i in range(1, 5)}
        assert len(h_ranks) == 1, "All H atoms in methane should be equivalent"


# ===================================================================
# H. Rotatable bonds
# ===================================================================


class TestRotatableBonds:
    def test_methane_no_rotatable(self, methane):
        assert methane.rotatable_bonds == []

    def test_ethanol_has_rotatable(self, ethanol):
        assert len(ethanol.rotatable_bonds) >= 1

    def test_rotatable_bonds_are_tuples(self, ethanol):
        for bond in ethanol.rotatable_bonds:
            assert len(bond) == 2

    def test_rotatable_bonds_zero_indexed(self, ethanol):
        for i, j in ethanol.rotatable_bonds:
            assert 0 <= i < ethanol.num_atoms
            assert 0 <= j < ethanol.num_atoms


# ===================================================================
# I. Immutable coordinate update
# ===================================================================


class TestCoordinateUpdate:
    def test_returns_new_molecule(self, methane):
        new_mol = methane.with_updated_coordinates(methane.coordinates)
        assert new_mol is not methane

    def test_coordinates_actually_updated(self, methane):
        coords = methane.coordinates.copy()
        coords[0] += 1.0
        new_mol = methane.with_updated_coordinates(coords)
        np.testing.assert_allclose(new_mol.coordinates[0], coords[0])

    def test_wrong_shape_raises(self, methane):
        with pytest.raises(ValueError, match="shape"):
            methane.with_updated_coordinates(np.zeros((3, 3)))


# ===================================================================
# J. Equality and repr
# ===================================================================


class TestMoleculeEquality:
    def test_same_smiles_equal(self):
        m1 = Molecule.from_smiles("CCO", embed_3d=False)
        m2 = Molecule.from_smiles("CCO", embed_3d=False)
        assert m1 == m2

    def test_different_smiles_not_equal(self):
        m1 = Molecule.from_smiles("CCO", embed_3d=False)
        m2 = Molecule.from_smiles("CCC", embed_3d=False)
        assert m1 != m2

    def test_not_equal_to_non_molecule(self, methane):
        assert methane != "not a molecule"


class TestMoleculeRepr:
    def test_repr_contains_name(self, methane):
        assert "methane" in repr(methane)

    def test_repr_contains_formula(self, methane):
        assert "CH4" in repr(methane)


# ===================================================================
# K. Config → Molecule integration
# ===================================================================


class TestConfigMoleculeIntegration:
    """Integration tests: config structure path loads a Molecule."""

    def test_config_structure_path_to_molecule(self, sdf_file):
        """A config pointing at an SDF should yield a valid Molecule."""
        cfg = PoltypeConfig(structure=str(sdf_file))
        mol = Molecule.from_file(cfg.structure)
        assert mol.num_atoms > 0
        assert mol.has_conformer

    def test_molecule_charge_matches_config(self, sdf_file):
        """Molecule formal charge should match config's total_charge default."""
        cfg = PoltypeConfig(structure=str(sdf_file), total_charge=0)
        mol = Molecule.from_file(cfg.structure)
        assert mol.formal_charge == cfg.total_charge

    def test_multiple_configs_independent(self):
        """Two PoltypeConfig instances should be independent."""
        c1 = PoltypeConfig(total_charge=0)
        c2 = PoltypeConfig(total_charge=-1)
        assert c1.total_charge != c2.total_charge

    def test_qm_config_overrides_propagate(self):
        """Sub-config overrides should be accessible through the parent."""
        cfg = PoltypeConfig(qm=QMConfig(opt_method="wB97X-D"))
        assert cfg.qm.opt_method == "wB97X-D"
        assert cfg.qm.opt_basis_set == "6-31G*"  # default preserved

    def test_resource_config_overrides_propagate(self):
        cfg = PoltypeConfig(resources=ResourceConfig(num_proc=16))
        assert cfg.resources.num_proc == 16
        assert cfg.resources.max_mem_gb == 16  # default preserved
