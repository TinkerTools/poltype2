"""
poltype.typing.type_record – per-atom type assignment record.

Each atom in a molecule receives a :class:`TypeRecord` that captures the
Tinker atom type, atom class, human-readable description, and the SMARTS
pattern that produced the match.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class TypeRecord:
    """Immutable record for a single atom's type assignment.

    Attributes
    ----------
    atom_index:
        0-based atom index (RDKit convention).
    atom_type:
        Tinker atom type integer.
    atom_class:
        Tinker atom class integer (may equal ``atom_type`` when no
        separate class mapping exists).
    description:
        Human-readable label, e.g. ``"Alkane CH3"``.
    smarts:
        The SMARTS pattern that matched this atom (empty string if the
        type was assigned by a fallback rule).
    """

    atom_index: int
    atom_type: int
    atom_class: int
    description: str = ""
    smarts: str = ""
