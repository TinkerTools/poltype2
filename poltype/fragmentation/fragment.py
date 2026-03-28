"""
poltype.fragmentation.fragment – data types for molecular fragments.

A :class:`Fragment` represents a chemically meaningful sub-structure
extracted from a parent molecule.  The :class:`FragmentResult`
collects all fragments produced by a single fragmentation run.

These dataclasses are plain containers that can be unit-tested in isolation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class Fragment:
    """A molecular fragment cut from a parent molecule.

    Attributes
    ----------
    fragment_id:
        Unique integer identifier within a fragmentation run
        (0-based).
    atom_map:
        Maps fragment atom index (0-based) → parent atom index
        (0-based).  This allows results computed on the fragment
        to be mapped back to the full molecule.
    rotatable_bond:
        The ``(i, j)`` pair of *parent* atom indices for the
        rotatable bond this fragment was generated to cover.
    cap_atoms:
        Set of fragment atom indices that are capping (hydrogen)
        atoms added during fragmentation (not present in the
        parent molecule).
    name:
        Human-readable label (e.g. ``"frag_0_bond_3_5"``).
    """

    fragment_id: int
    atom_map: Dict[int, int]
    rotatable_bond: Tuple[int, int]
    cap_atoms: List[int] = field(default_factory=list)
    name: str = ""

    @property
    def num_atoms(self) -> int:
        """Number of atoms in this fragment (including caps)."""
        return len(self.atom_map)

    @property
    def parent_atom_indices(self) -> List[int]:
        """Sorted list of parent atom indices covered by this fragment."""
        return sorted(self.atom_map.values())

    @property
    def num_cap_atoms(self) -> int:
        """Number of capping atoms added during fragmentation."""
        return len(self.cap_atoms)


@dataclass
class FragmentResult:
    """Collection of fragments produced by a single fragmentation run.

    Attributes
    ----------
    fragments:
        Ordered list of :class:`Fragment` objects.
    parent_molecule_name:
        Name of the parent molecule that was fragmented.
    fragmentation_method:
        Name of the method used (e.g. ``"wbo"``).
    metadata:
        Arbitrary extra information from the fragmentation
        (e.g. WBO threshold, growth cycles).
    """

    fragments: List[Fragment] = field(default_factory=list)
    parent_molecule_name: str = ""
    fragmentation_method: str = ""
    metadata: Dict[str, object] = field(default_factory=dict)

    @property
    def num_fragments(self) -> int:
        """Number of fragments in this result."""
        return len(self.fragments)

    @property
    def total_fragment_atoms(self) -> int:
        """Total number of atoms across all fragments."""
        return sum(f.num_atoms for f in self.fragments)

    @property
    def covered_parent_atoms(self) -> set:
        """Set of parent atom indices covered by at least one fragment."""
        covered: set = set()
        for frag in self.fragments:
            covered.update(frag.parent_atom_indices)
        return covered
