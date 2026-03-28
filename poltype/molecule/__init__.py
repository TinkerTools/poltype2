"""
poltype.molecule – canonical molecule representation.

Re-exports the :class:`Molecule` wrapper so callers can write::

    from poltype.molecule import Molecule
"""

from __future__ import annotations

from poltype.molecule.molecule import Molecule

__all__ = ["Molecule"]
