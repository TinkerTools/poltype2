"""
tests/unit/test_molecule.py – unit tests for poltype.molecule.Molecule.

Tests cover:
- Construction from SMILES.
- Basic property accessors (num_atoms, formula, formal_charge, etc.).
- Index conversion helpers.
- Symmetry class computation.
- Rotatable bond detection.
- Coordinate update returning a new Molecule.
- Equality comparison.
- Repr.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from poltype.molecule.molecule import Molecule


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def methane() -> Molecule:
    """CH4 – 5 atoms, no rotatable bonds, single symmetry class for H."""
    return Molecule.from_smiles("C", name="methane", embed_3d=True)


@pytest.fixture
def ethanol() -> Molecule:
    """Ethanol (with explicit H) – has one C-C rotatable bond."""
    return Molecule.from_smiles("CCO", name="ethanol", embed_3d=True)


@pytest.fixture
def aspirin_sdf(tmp_path) -> Path:
    """Write a valid aspirin SDF (RDKit-generated) to a temp file."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    sdf_block = Chem.MolToMolBlock(mol)
    p = tmp_path / "aspirin.sdf"
    p.write_text(sdf_block + "\n$$$$\n")
    return p


# ---------------------------------------------------------------------------
# Construction tests
# ---------------------------------------------------------------------------


class TestMoleculeConstruction:
    def test_from_smiles_basic(self):
        mol = Molecule.from_smiles("C", embed_3d=False)
        assert mol is not None

    def test_from_smiles_name(self):
        mol = Molecule.from_smiles("C", name="methane", embed_3d=False)
        assert mol.name == "methane"

    def test_from_smiles_with_embed(self, methane):
        assert methane.has_conformer

    def test_from_rdmol_direct(self):
        rdmol = Chem.MolFromSmiles("C")
        rdmol = Chem.AddHs(rdmol)
        mol = Molecule(rdmol, name="methane_direct")
        assert mol.num_atoms == 5

    def test_from_file_sdf(self, aspirin_sdf):
        mol = Molecule.from_file(aspirin_sdf)
        assert mol is not None
        assert mol.num_atoms > 0

    def test_from_file_default_name(self, aspirin_sdf):
        mol = Molecule.from_file(aspirin_sdf)
        assert mol.name == "aspirin"

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


# ---------------------------------------------------------------------------
# Property tests
# ---------------------------------------------------------------------------


class TestMoleculeProperties:
    def test_num_atoms_methane(self, methane):
        # CH4 with explicit H = 5
        assert methane.num_atoms == 5

    def test_num_heavy_atoms_methane(self, methane):
        assert methane.num_heavy_atoms == 1

    def test_formal_charge_neutral(self, methane):
        assert methane.formal_charge == 0

    def test_formal_charge_ion(self):
        mol = Molecule.from_smiles("[NH4+]", embed_3d=False)
        assert mol.formal_charge == 1

    def test_formula_methane(self, methane):
        assert methane.formula == "CH4"

    def test_has_conformer_after_embed(self, methane):
        assert methane.has_conformer is True

    def test_has_conformer_no_embed(self):
        mol = Molecule.from_smiles("C", embed_3d=False)
        assert mol.has_conformer is False

    def test_coordinates_shape(self, methane):
        coords = methane.coordinates
        assert coords.shape == (5, 3)

    def test_coordinates_no_conformer_raises(self):
        mol = Molecule.from_smiles("C", embed_3d=False)
        with pytest.raises(RuntimeError, match="no conformer"):
            _ = mol.coordinates

    def test_rdmol_is_chem_mol(self, methane):
        assert isinstance(methane.rdmol, Chem.Mol)

    def test_work_dir_is_absolute(self, methane):
        assert methane.work_dir.is_absolute()


# ---------------------------------------------------------------------------
# Index conversion tests
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Symmetry class tests
# ---------------------------------------------------------------------------


class TestSymmetryClasses:
    def test_returns_dict(self, methane):
        sym = methane.symmetry_classes
        assert isinstance(sym, dict)

    def test_all_atoms_present(self, methane):
        sym = methane.symmetry_classes
        assert set(sym.keys()) == set(range(methane.num_atoms))

    def test_equivalent_hydrogens(self, methane):
        # All 4 H atoms in methane should have the same symmetry rank.
        sym = methane.symmetry_classes
        # Atom 0 is C; atoms 1-4 are H
        h_ranks = {sym[i] for i in range(1, 5)}
        assert len(h_ranks) == 1, "All H atoms should be equivalent"


# ---------------------------------------------------------------------------
# Rotatable bond tests
# ---------------------------------------------------------------------------


class TestRotatableBonds:
    def test_methane_no_rotatable_bonds(self, methane):
        assert methane.rotatable_bonds == []

    def test_ethanol_has_rotatable_bond(self, ethanol):
        # C-C and C-O bonds in CCO should be rotatable
        assert len(ethanol.rotatable_bonds) >= 1

    def test_rotatable_bonds_are_tuples(self, ethanol):
        for bond in ethanol.rotatable_bonds:
            assert len(bond) == 2

    def test_rotatable_bonds_zero_indexed(self, ethanol):
        for i, j in ethanol.rotatable_bonds:
            assert i >= 0 and j >= 0
            assert i < ethanol.num_atoms and j < ethanol.num_atoms


# ---------------------------------------------------------------------------
# Coordinate update tests
# ---------------------------------------------------------------------------


class TestWithUpdatedCoordinates:
    def test_returns_new_molecule(self, methane):
        coords = methane.coordinates
        new_mol = methane.with_updated_coordinates(coords)
        assert new_mol is not methane

    def test_coordinates_updated(self, methane):
        coords = methane.coordinates.copy()
        coords[0] += 1.0  # shift first atom
        new_mol = methane.with_updated_coordinates(coords)
        np.testing.assert_allclose(new_mol.coordinates[0], coords[0])

    def test_wrong_shape_raises(self, methane):
        bad_coords = np.zeros((3, 3))
        with pytest.raises(ValueError, match="shape"):
            methane.with_updated_coordinates(bad_coords)


# ---------------------------------------------------------------------------
# Equality and repr tests
# ---------------------------------------------------------------------------


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
        r = repr(methane)
        assert "methane" in r

    def test_repr_contains_formula(self, methane):
        r = repr(methane)
        assert "CH4" in r
