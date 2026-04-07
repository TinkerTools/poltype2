"""Shared helpers for materialized fitted parameter artifacts."""

from __future__ import annotations


def atom_type_for_index(atom_index: int, type_map: dict[int, int]) -> int:
    """Return a stable atom type for an atom index."""
    atom_type = type_map.get(atom_index, atom_index + 1)
    return atom_type if atom_type != 0 else atom_index + 1


def torsion_indices_for_bond(
    rdmol,
    bond: tuple[int, int],
) -> tuple[int, int, int, int] | None:
    """Build a simple four-atom torsion tuple for a central bond."""
    j, k = bond
    left = sorted(
        neighbor.GetIdx()
        for neighbor in rdmol.GetAtomWithIdx(j).GetNeighbors()
        if neighbor.GetIdx() != k
    )
    right = sorted(
        neighbor.GetIdx()
        for neighbor in rdmol.GetAtomWithIdx(k).GetNeighbors()
        if neighbor.GetIdx() != j
    )
    if not left or not right:
        return None
    return (left[0], j, k, right[0])
