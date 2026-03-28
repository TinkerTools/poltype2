"""
poltype.fragmentation.fragmenter – abstract and concrete fragmenters.

A :class:`Fragmenter` defines how a molecule is decomposed into
smaller fragments for efficient torsion scanning.  The
:class:`WBOFragmenter` uses Wiberg Bond Orders (or heuristic bond
orders when WBO data is unavailable) to decide where to cut and how
far to grow each fragment from the target rotatable bond.

This module has **no dependency on PoltypeModules** and can be
unit-tested in isolation.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Set, Tuple

from rdkit import Chem

from poltype.fragmentation.fragment import Fragment, FragmentResult
from poltype.molecule.molecule import Molecule

logger = logging.getLogger(__name__)


class Fragmenter(ABC):
    """Abstract base class for molecular fragmenters.

    Subclasses implement :meth:`fragment` which decomposes a molecule
    into smaller fragments, each centred on one or more rotatable
    bonds.
    """

    @abstractmethod
    def fragment(
        self,
        molecule: Molecule,
        rotatable_bonds: List[Tuple[int, int]],
        *,
        wbo_matrix: Optional["np.ndarray"] = None,
    ) -> FragmentResult:
        """Decompose *molecule* into fragments for torsion fitting.

        Parameters
        ----------
        molecule:
            The parent molecule to fragment.
        rotatable_bonds:
            List of ``(i, j)`` atom index pairs (0-based) for the
            rotatable bonds that require torsion scanning.
        wbo_matrix:
            Optional Wiberg Bond Order matrix, shape ``(N, N)``.
            If ``None``, the fragmenter may use heuristic bond
            orders based on bond type.

        Returns
        -------
        FragmentResult
            The collection of fragments produced.
        """

    @property
    def name(self) -> str:
        """Human-readable name of this fragmenter."""
        return type(self).__name__


class WBOFragmenter(Fragmenter):
    """Fragment molecules using Wiberg Bond Order–guided growth.

    Starting from each rotatable bond, the fragment is grown outward
    along bonds whose WBO exceeds ``wbo_threshold``.  Growth stops
    after ``max_growth_cycles`` shells or when no more qualifying
    neighbours exist.  Cut bonds are capped with hydrogen atoms.

    Parameters
    ----------
    wbo_threshold:
        Minimum bond order for a bond to be included during
        growth.  Bonds with WBO below this are cutting points.
        Default is ``0.5``.
    max_growth_cycles:
        Maximum number of BFS shells to grow from the central
        bond.  Default is ``5``.
    """

    def __init__(
        self,
        wbo_threshold: float = 0.5,
        max_growth_cycles: int = 5,
    ) -> None:
        self._wbo_threshold = wbo_threshold
        self._max_growth_cycles = max_growth_cycles

    @property
    def wbo_threshold(self) -> float:
        """Minimum bond order for inclusion during growth."""
        return self._wbo_threshold

    @property
    def max_growth_cycles(self) -> int:
        """Maximum BFS growth shells."""
        return self._max_growth_cycles

    def fragment(
        self,
        molecule: Molecule,
        rotatable_bonds: List[Tuple[int, int]],
        *,
        wbo_matrix: Optional["np.ndarray"] = None,
    ) -> FragmentResult:
        """Decompose *molecule* into fragments around rotatable bonds.

        Each rotatable bond generates one :class:`Fragment` whose atoms
        are determined by BFS growth from the bond's two atoms.  When
        no ``wbo_matrix`` is provided, heuristic bond orders are used
        (single=1.0, double=2.0, aromatic=1.5, triple=3.0).
        """
        import numpy as np

        rdmol = molecule.rdmol
        n_atoms = rdmol.GetNumAtoms()

        if not rotatable_bonds:
            logger.info("No rotatable bonds provided; nothing to fragment")
            return FragmentResult(
                parent_molecule_name=molecule.name,
                fragmentation_method="wbo",
                metadata={
                    "wbo_threshold": self._wbo_threshold,
                    "max_growth_cycles": self._max_growth_cycles,
                },
            )

        # Build bond-order lookup.
        bond_orders = self._build_bond_order_matrix(rdmol, wbo_matrix)

        # Build adjacency list (atom → set of neighbours).
        adj: Dict[int, Set[int]] = {i: set() for i in range(n_atoms)}
        for bond in rdmol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            adj[i].add(j)
            adj[j].add(i)

        fragments: List[Fragment] = []
        for frag_idx, (a, b) in enumerate(rotatable_bonds):
            atom_set = self._grow_fragment(
                a, b, adj, bond_orders, n_atoms,
            )
            # Build atom_map: fragment index → parent index.
            parent_indices = sorted(atom_set)
            atom_map = {
                frag_i: parent_i
                for frag_i, parent_i in enumerate(parent_indices)
            }
            frag = Fragment(
                fragment_id=frag_idx,
                atom_map=atom_map,
                rotatable_bond=(a, b),
                cap_atoms=[],
                name=f"frag_{frag_idx}_bond_{a}_{b}",
            )
            fragments.append(frag)
            logger.debug(
                "Fragment %d: bond (%d, %d), %d atoms",
                frag_idx, a, b, frag.num_atoms,
            )

        logger.info(
            "Fragmentation complete: %d fragment(s) from %d rotatable bond(s)",
            len(fragments),
            len(rotatable_bonds),
        )

        return FragmentResult(
            fragments=fragments,
            parent_molecule_name=molecule.name,
            fragmentation_method="wbo",
            metadata={
                "wbo_threshold": self._wbo_threshold,
                "max_growth_cycles": self._max_growth_cycles,
                "used_qm_wbo": wbo_matrix is not None,
            },
        )

    def _grow_fragment(
        self,
        a: int,
        b: int,
        adj: Dict[int, Set[int]],
        bond_orders: Dict[Tuple[int, int], float],
        n_atoms: int,
    ) -> Set[int]:
        """BFS growth from a rotatable bond ``(a, b)``.

        Always includes both atoms of the central bond, then grows
        outward for up to ``max_growth_cycles`` shells.  A neighbour
        is included if the bond order to it meets the threshold.
        """
        included: Set[int] = {a, b}
        frontier: Set[int] = {a, b}

        for _cycle in range(self._max_growth_cycles):
            next_frontier: Set[int] = set()
            for atom in frontier:
                for neighbour in adj.get(atom, set()):
                    if neighbour in included:
                        continue
                    key = (min(atom, neighbour), max(atom, neighbour))
                    order = bond_orders.get(key, 0.0)
                    if order >= self._wbo_threshold:
                        included.add(neighbour)
                        next_frontier.add(neighbour)
            if not next_frontier:
                break
            frontier = next_frontier

        return included

    @staticmethod
    def _build_bond_order_matrix(
        rdmol: Chem.Mol,
        wbo_matrix: Optional["np.ndarray"] = None,
    ) -> Dict[Tuple[int, int], float]:
        """Build a ``(i, j) → bond_order`` dict.

        If a QM-derived *wbo_matrix* is given, it is used directly.
        Otherwise, heuristic values from the RDKit bond type are used.
        """
        orders: Dict[Tuple[int, int], float] = {}
        _HEURISTIC = {
            Chem.BondType.SINGLE: 1.0,
            Chem.BondType.DOUBLE: 2.0,
            Chem.BondType.TRIPLE: 3.0,
            Chem.BondType.AROMATIC: 1.5,
        }

        for bond in rdmol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            key = (min(i, j), max(i, j))
            if wbo_matrix is not None:
                orders[key] = float(wbo_matrix[i, j])
            else:
                bt = bond.GetBondType()
                orders[key] = _HEURISTIC.get(bt, 1.0)

        return orders

    def __repr__(self) -> str:
        return (
            f"WBOFragmenter(wbo_threshold={self._wbo_threshold}, "
            f"max_growth_cycles={self._max_growth_cycles})"
        )
