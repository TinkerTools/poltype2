"""
poltype.database.match_result – result of a parameter database lookup.

A :class:`DatabaseMatchResult` summarises which atoms (or terms) in a
molecule already have parameters in the library, so the pipeline can
decide which stages to run and which to skip.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set


@dataclass
class DatabaseMatchResult:
    """Result of matching a molecule against a parameter database.

    Attributes
    ----------
    matched_atom_indices:
        0-based atom indices for which types were found in the database.
    unmatched_atom_indices:
        0-based atom indices with no database match.
    matched_torsions:
        Set of 4-tuples ``(i, j, k, l)`` (0-based) for torsion terms
        already in the database.
    source:
        Name of the database that produced this result.
    metadata:
        Arbitrary extra information from the lookup (e.g. database
        version, number of rules tried).
    """

    matched_atom_indices: List[int] = field(default_factory=list)
    unmatched_atom_indices: List[int] = field(default_factory=list)
    matched_torsions: Set[tuple] = field(default_factory=set)
    source: str = ""
    metadata: Dict[str, object] = field(default_factory=dict)

    @property
    def all_atoms_matched(self) -> bool:
        """``True`` if every atom was found in the database."""
        return len(self.unmatched_atom_indices) == 0

    @property
    def coverage(self) -> float:
        """Fraction of atoms matched (0.0 – 1.0)."""
        total = len(self.matched_atom_indices) + len(self.unmatched_atom_indices)
        if total == 0:
            return 0.0
        return len(self.matched_atom_indices) / total
