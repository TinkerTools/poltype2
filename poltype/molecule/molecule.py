"""
poltype.molecule.molecule – the canonical ``Molecule`` wrapper.

Design goals
------------
* **One representation** – RDKit is the primary in-memory format.
* **Lazy OpenBabel** – ``ob_mol`` is computed on first access and cached;
  it is an *implementation detail* for code that still requires OB.
* **0-based indices everywhere** – explicit conversion helpers remove the
  silent ±1 arithmetic scattered throughout the legacy code.
* **Absolute paths** – every Molecule knows its working directory; no
  ``os.chdir`` is needed.
* **Immutable updates** – methods that modify geometry return a *new*
  ``Molecule`` rather than mutating in place.

This module has **no dependency on PoltypeModules** and can be unit-tested
in isolation.

Example::

    from poltype.molecule import Molecule
    mol = Molecule.from_file("aspirin.sdf")
    print(mol.num_atoms)           # 21 (RDKit, 0-based)
    print(mol.to_tinker_idx(0))    # 1  (Tinker, 1-based)
    sym = mol.symmetry_classes     # {0: 1, 1: 2, …}
"""

from __future__ import annotations

import copy
from functools import cached_property
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# RDKit is a hard dependency for the new architecture.
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, rdmolops


class Molecule:
    """Canonical molecule representation for the Poltype2 pipeline.

    Parameters
    ----------
    rdmol:
        RDKit molecule (should have a conformer if 3-D coordinates are
        needed).
    work_dir:
        Absolute path to the working directory for this molecule.  All
        file I/O performed on behalf of this molecule should use paths
        relative to ``work_dir``.
    name:
        Optional human-readable name (used in log messages).
    """

    def __init__(
        self,
        rdmol: Chem.Mol,
        work_dir: Path = Path("."),
        name: str = "",
    ) -> None:
        if rdmol is None:
            raise ValueError("rdmol must not be None")
        # Store a read-only copy so callers cannot mutate our internal state.
        self._rdmol: Chem.Mol = copy.deepcopy(rdmol)
        self.work_dir: Path = Path(work_dir).resolve()
        self.name: str = name

    # ------------------------------------------------------------------
    # Alternate constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        work_dir: Optional[Path] = None,
        name: str = "",
    ) -> "Molecule":
        """Create a ``Molecule`` from an SDF, MOL, or MOL2 file.

        Parameters
        ----------
        path:
            Path to the input molecule file.
        work_dir:
            Working directory.  Defaults to the directory containing
            *path*.
        name:
            Optional molecule name override.  If empty, the file stem is
            used.
        """
        path = Path(path).resolve()
        if not path.exists():
            raise FileNotFoundError(f"Molecule file not found: {path}")

        suffix = path.suffix.lower()
        rdmol: Optional[Chem.Mol] = None

        if suffix in (".sdf", ".mol"):
            supplier = Chem.SDMolSupplier(str(path), removeHs=False)
            for mol in supplier:
                if mol is not None:
                    rdmol = mol
                    break
        elif suffix == ".mol2":
            rdmol = Chem.MolFromMol2File(str(path), removeHs=False)
        else:
            raise ValueError(
                f"Unsupported file format '{suffix}'. "
                "Supported: .sdf, .mol, .mol2"
            )

        if rdmol is None:
            raise ValueError(f"Failed to parse molecule from {path}")

        mol_name = name or path.stem
        wd = work_dir if work_dir is not None else path.parent
        return cls(rdmol, work_dir=wd, name=mol_name)

    @classmethod
    def from_smiles(
        cls,
        smiles: str,
        work_dir: Path = Path("."),
        name: str = "",
        add_hydrogens: bool = True,
        embed_3d: bool = True,
    ) -> "Molecule":
        """Create a ``Molecule`` from a SMILES string.

        Parameters
        ----------
        smiles:
            SMILES string for the molecule.
        work_dir:
            Working directory.
        name:
            Optional molecule name.
        add_hydrogens:
            If ``True`` (default), explicit hydrogens are added before
            embedding.
        embed_3d:
            If ``True`` (default), an initial 3-D conformer is generated
            using RDKit's ETKDG.
        """
        rdmol = Chem.MolFromSmiles(smiles)
        if rdmol is None:
            raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
        if add_hydrogens:
            rdmol = Chem.AddHs(rdmol)
        if embed_3d:
            AllChem.EmbedMolecule(rdmol, AllChem.ETKDGv3())
        return cls(rdmol, work_dir=work_dir, name=name)

    # ------------------------------------------------------------------
    # Core properties
    # ------------------------------------------------------------------

    @property
    def rdmol(self) -> Chem.Mol:
        """Read-only access to the underlying RDKit molecule."""
        return self._rdmol

    @property
    def num_atoms(self) -> int:
        """Total number of atoms (including hydrogens)."""
        return self._rdmol.GetNumAtoms()

    @property
    def num_heavy_atoms(self) -> int:
        """Number of heavy (non-hydrogen) atoms."""
        return rdMolDescriptors.CalcNumHeavyAtoms(self._rdmol)

    @property
    def formal_charge(self) -> int:
        """Net formal charge of the molecule."""
        return Chem.GetFormalCharge(self._rdmol)

    @property
    def formula(self) -> str:
        """Molecular formula string (e.g. ``"C9H8O4"``)."""
        return rdMolDescriptors.CalcMolFormula(self._rdmol)

    @property
    def has_conformer(self) -> bool:
        """Return ``True`` if at least one 3-D conformer is present."""
        return self._rdmol.GetNumConformers() > 0

    @property
    def coordinates(self) -> np.ndarray:
        """Atomic coordinates of the first conformer, shape ``(N, 3)``."""
        if not self.has_conformer:
            raise RuntimeError(
                "Molecule has no conformer; embed or load a 3-D structure first."
            )
        conf = self._rdmol.GetConformer(0)
        return np.array(
            [conf.GetAtomPosition(i) for i in range(self.num_atoms)]
        )

    # ------------------------------------------------------------------
    # Lazy / cached properties
    # ------------------------------------------------------------------

    @cached_property
    def symmetry_classes(self) -> Dict[int, int]:
        """Map from 0-based atom index to canonical symmetry class integer.

        Uses RDKit's Morgan algorithm (radius 0 = topological equivalence).
        Two atoms in the same symmetry class have the same chemical
        environment and can be treated as equivalent for parameter
        assignment.
        """
        classes = list(
            rdmolops.GetSymmSSSR(self._rdmol)  # ensures ring info is cached
        )
        ranks = list(
            Chem.CanonicalRankAtoms(self._rdmol, breakTies=False)
        )
        return {i: int(r) for i, r in enumerate(ranks)}

    @cached_property
    def rotatable_bonds(self) -> List[Tuple[int, int]]:
        """List of ``(i, j)`` pairs (0-based) for rotatable bonds.

        Rotatable bonds are those that are:
        * Single bonds
        * Not in a ring
        * Not terminal (both atoms have degree > 1)
        """
        bonds = []
        for bond in self._rdmol.GetBonds():
            if bond.GetBondTypeAsDouble() != 1.0:
                continue
            if bond.IsInRing():
                continue
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            ai = self._rdmol.GetAtomWithIdx(i)
            aj = self._rdmol.GetAtomWithIdx(j)
            if ai.GetDegree() > 1 and aj.GetDegree() > 1:
                bonds.append((i, j))
        return bonds

    @cached_property
    def ob_mol(self) -> "openbabel.OBMol":  # type: ignore[name-defined]
        """Lazy OpenBabel representation (computed on first access).

        This property exists purely as a compatibility bridge for legacy
        code that requires an OpenBabel ``OBMol``.  New code should work
        exclusively with :attr:`rdmol`.

        Raises
        ------
        ImportError
            If OpenBabel is not installed in the current environment.
        RuntimeError
            If the conversion from RDKit to OpenBabel fails.
        """
        try:
            from openbabel import openbabel as ob
            from rdkit.Chem import AllChem
        except ImportError as exc:  # pragma: no cover
            raise ImportError(
                "OpenBabel is required for ob_mol conversion."
            ) from exc

        # Write to a temporary in-memory SDF string and re-read with OB.
        sdf_block = Chem.MolToMolBlock(self._rdmol)
        if sdf_block is None:
            raise RuntimeError(
                "RDKit could not serialise molecule to molblock."
            )
        obconv = ob.OBConversion()
        obconv.SetInAndOutFormats("mol", "mol")
        ob_mol = ob.OBMol()
        if not obconv.ReadString(ob_mol, sdf_block):
            raise RuntimeError(
                "OpenBabel could not parse the RDKit-generated molblock."
            )
        return ob_mol

    # ------------------------------------------------------------------
    # Index conversion helpers
    # ------------------------------------------------------------------

    @staticmethod
    def to_tinker_idx(rdkit_idx: int) -> int:
        """Convert a 0-based RDKit atom index to a 1-based Tinker index."""
        return rdkit_idx + 1

    @staticmethod
    def to_rdkit_idx(tinker_idx: int) -> int:
        """Convert a 1-based Tinker atom index to a 0-based RDKit index."""
        return tinker_idx - 1

    # ------------------------------------------------------------------
    # Geometry update (returns a *new* Molecule)
    # ------------------------------------------------------------------

    def with_updated_coordinates(self, coords: np.ndarray) -> "Molecule":
        """Return a new ``Molecule`` with the given atomic coordinates.

        Parameters
        ----------
        coords:
            New coordinates array, shape ``(N, 3)``, where N is
            :attr:`num_atoms`.

        Returns
        -------
        Molecule
            A new ``Molecule`` instance with the same topology but
            updated 3-D coordinates.
        """
        if coords.shape != (self.num_atoms, 3):
            raise ValueError(
                f"coords shape {coords.shape} does not match "
                f"num_atoms={self.num_atoms}"
            )
        new_rdmol = copy.deepcopy(self._rdmol)
        if new_rdmol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(new_rdmol, AllChem.ETKDGv3())
        conf = new_rdmol.GetConformer(0)
        from rdkit.Geometry import Point3D
        for i, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
        return Molecule(new_rdmol, work_dir=self.work_dir, name=self.name)

    # ------------------------------------------------------------------
    # Dunder methods
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"Molecule(name={self.name!r}, formula={self.formula}, "
            f"atoms={self.num_atoms}, work_dir={self.work_dir})"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Molecule):
            return NotImplemented
        # Compare by canonical SMILES (ignores conformer differences).
        return (
            Chem.MolToSmiles(self._rdmol)
            == Chem.MolToSmiles(other._rdmol)
        )
