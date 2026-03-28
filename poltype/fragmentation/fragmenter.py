"""
poltype.fragmentation.fragmenter – abstract and concrete fragmenters.

A :class:`Fragmenter` defines how a molecule is decomposed into
smaller fragments for efficient torsion scanning.

* :class:`WBOFragmenter` uses Wiberg Bond Orders (or heuristic bond
  orders when WBO data is unavailable) to decide where to cut and
  how far to grow each fragment from the target rotatable bond.
* :class:`RuleBasedFragmenter` uses predefined chemical rules
  (similar to the post-processing script
  ``PoltypeModules/lTorsionFragmentPostProcessing.py``) to trim
  fragments without requiring WBO data.

This module has **no dependency on PoltypeModules** and can be
unit-tested in isolation.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Dict, List, Optional, Set, Tuple

from rdkit import Chem
from rdkit.Chem import RingInfo

from poltype.fragmentation.fragment import Fragment, FragmentResult
from poltype.molecule.molecule import Molecule

if TYPE_CHECKING:
    import numpy as np

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
        wbo_matrix: Optional[np.ndarray] = None,
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
        wbo_matrix: Optional[np.ndarray] = None,
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
        wbo_matrix: Optional[np.ndarray] = None,
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


class RuleBasedFragmenter(Fragmenter):
    """Fragment molecules using predefined chemical rules.

    Instead of relying on Wiberg Bond Order (WBO) data, this
    fragmenter applies explicit chemical rules inspired by the
    post-processing script ``lTorsionFragmentPostProcessing.py``
    to trim fragments around each rotatable bond.

    Key rules
    ---------
    1. **Eligible bond cutting** – all single bonds between heavy
       atoms and C-C double/aromatic bonds may be cut.
    2. **Nitrogen protection** – bonds involving divalent nitrogen
       (``[#7X2]``) are never cut.
    3. **Torsion region protection** – bonds whose *both* atoms
       belong to the immediate neighbourhood of the rotatable bond
       (the two central atoms and all their direct neighbours) are
       never cut.
    4. **Aromatic substituent removal** – non-ring neighbours of
       aromatic ring atoms that are outside the torsion region are
       removed.
    5. **Fused-ring trimming** – distant fused rings (connected
       via C-C bonds in fused positions) are removed when they
       do not share a ring with the torsion atoms.
    6. **Connectivity cleanup** – after all cuts, connected
       components that do not contain the torsion region are
       discarded.

    Parameters
    ----------
    max_fragment_atoms:
        Soft upper bound on fragment size.  If a fragment exceeds
        this count, extra trimming heuristics may be applied in
        future versions.  Currently informational only.  Default
        is ``50``.
    """

    def __init__(self, max_fragment_atoms: int = 50) -> None:
        self._max_fragment_atoms = max_fragment_atoms

    @property
    def max_fragment_atoms(self) -> int:
        """Soft upper bound on fragment atom count."""
        return self._max_fragment_atoms

    # -----------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------

    def fragment(
        self,
        molecule: Molecule,
        rotatable_bonds: List[Tuple[int, int]],
        *,
        wbo_matrix: Optional["np.ndarray"] = None,
    ) -> FragmentResult:
        """Decompose *molecule* into fragments using chemical rules.

        The ``wbo_matrix`` parameter is accepted for interface
        compatibility but is **ignored** – this fragmenter does not
        use WBO data.
        """
        rdmol = molecule.rdmol

        if not rotatable_bonds:
            logger.info("No rotatable bonds provided; nothing to fragment")
            return FragmentResult(
                parent_molecule_name=molecule.name,
                fragmentation_method="rule_based",
                metadata={"max_fragment_atoms": self._max_fragment_atoms},
            )

        fragments: List[Fragment] = []
        for frag_idx, (a, b) in enumerate(rotatable_bonds):
            atom_set = self._rule_based_fragment(rdmol, a, b)
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
            "Rule-based fragmentation complete: %d fragment(s) "
            "from %d rotatable bond(s)",
            len(fragments),
            len(rotatable_bonds),
        )
        return FragmentResult(
            fragments=fragments,
            parent_molecule_name=molecule.name,
            fragmentation_method="rule_based",
            metadata={"max_fragment_atoms": self._max_fragment_atoms},
        )

    # -----------------------------------------------------------------
    # Private helpers
    # -----------------------------------------------------------------

    def _rule_based_fragment(
        self,
        rdmol: Chem.Mol,
        a: int,
        b: int,
    ) -> Set[int]:
        """Determine the set of parent atom indices to keep for bond ``(a, b)``.

        Applies the rule-based trimming algorithm.
        """
        n_atoms = rdmol.GetNumAtoms()
        all_atoms = set(range(n_atoms))

        # 1. Build torsion region: atoms a, b and all their direct neighbours.
        torsion_region = self._torsion_region(rdmol, a, b)

        # 2. Identify eligible bonds (all heavy–heavy single bonds + C~C bonds).
        eligible_bonds = self._eligible_bonds(rdmol)

        # 3. Identify protected bonds (divalent nitrogen).
        special_nitrogen_bonds = self._special_nitrogen_bonds(rdmol)

        # 4. Build adjacency as a mutable edge set for cutting.
        edges: Set[Tuple[int, int]] = set()
        adj: Dict[int, Set[int]] = {i: set() for i in range(n_atoms)}
        for bond in rdmol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            key = (min(i, j), max(i, j))
            edges.add(key)
            adj[i].add(j)
            adj[j].add(i)

        # 5. Identify substituent bonds to cut.
        bonds_to_cut: List[Tuple[int, int]] = []
        ring_info = rdmol.GetRingInfo()

        for idx in range(n_atoms):
            atom = rdmol.GetAtomWithIdx(idx)
            in_ring = atom.IsInRing()
            is_aromatic = atom.GetIsAromatic()

            if in_ring and is_aromatic:
                for neig in atom.GetNeighbors():
                    neig_idx = neig.GetIdx()
                    same_ring = RingInfo.AreAtomsInSameRing(
                        ring_info, idx, neig_idx,
                    )
                    if same_ring or neig.GetAtomicNum() == 1:
                        continue
                    pair = (min(idx, neig_idx), max(idx, neig_idx))
                    if (
                        pair in eligible_bonds
                        and pair not in special_nitrogen_bonds
                        and not {pair[0], pair[1]}.issubset(torsion_region)
                        and pair not in bonds_to_cut
                    ):
                        bonds_to_cut.append(pair)
            else:
                for neig in atom.GetNeighbors():
                    neig_idx = neig.GetIdx()
                    pair = (min(idx, neig_idx), max(idx, neig_idx))
                    if (
                        pair in eligible_bonds
                        and pair not in special_nitrogen_bonds
                        and not {pair[0], pair[1]}.issubset(torsion_region)
                        and pair not in bonds_to_cut
                    ):
                        bonds_to_cut.append(pair)

        # 6. Cut the identified bonds from adjacency.
        for i, j in bonds_to_cut:
            adj[i].discard(j)
            adj[j].discard(i)

        # 7. Find connected components; keep only the one containing
        #    the torsion region.
        atoms_to_remove: Set[int] = set()
        components = self._connected_components(adj, n_atoms)
        for comp in components:
            if not torsion_region.issubset(comp):
                atoms_to_remove.update(comp)

        # 8. Handle fused rings: remove atoms on distant fused rings.
        atoms_to_remove.update(
            self._fused_ring_atoms_to_remove(rdmol, torsion_region),
        )

        # 9. Remove isolated atoms from adjacency and re-check components.
        for idx in atoms_to_remove:
            for neig in list(adj.get(idx, set())):
                adj[neig].discard(idx)
            adj[idx] = set()
        components = self._connected_components(adj, n_atoms)
        for comp in components:
            if not torsion_region.issubset(comp):
                atoms_to_remove.update(comp)

        # 10. Safety: never remove torsion region atoms.
        if atoms_to_remove & torsion_region:
            atoms_to_remove.clear()

        kept = all_atoms - atoms_to_remove
        return kept

    @staticmethod
    def _torsion_region(rdmol: Chem.Mol, a: int, b: int) -> Set[int]:
        """Return atoms ``a``, ``b`` and all their direct heavy-atom neighbours."""
        region: Set[int] = {a, b}
        for center in (a, b):
            atom = rdmol.GetAtomWithIdx(center)
            for neig in atom.GetNeighbors():
                region.add(neig.GetIdx())
        return region

    @staticmethod
    def _eligible_bonds(rdmol: Chem.Mol) -> Set[Tuple[int, int]]:
        """Bonds eligible for cutting: heavy–heavy single + C~C bonds."""
        eligible: Set[Tuple[int, int]] = set()

        # All bonds between non-hydrogen atoms.
        for bond in rdmol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            if (
                rdmol.GetAtomWithIdx(i).GetAtomicNum() != 1
                and rdmol.GetAtomWithIdx(j).GetAtomicNum() != 1
            ):
                eligible.add((min(i, j), max(i, j)))

        # C~C double/aromatic bonds (already captured above if
        # both are non-H, but this ensures all C-C bond types).
        for bond in rdmol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            if (
                rdmol.GetAtomWithIdx(i).GetAtomicNum() == 6
                and rdmol.GetAtomWithIdx(j).GetAtomicNum() == 6
            ):
                eligible.add((min(i, j), max(i, j)))

        return eligible

    @staticmethod
    def _special_nitrogen_bonds(rdmol: Chem.Mol) -> Set[Tuple[int, int]]:
        """Bonds involving divalent nitrogen that must not be cut."""
        protected: Set[Tuple[int, int]] = set()

        # [#7X2][*] – nitrogen with exactly 2 connections.
        pattern = Chem.MolFromSmarts("[#7X2][*]")
        if pattern is not None:
            for match in rdmol.GetSubstructMatches(pattern, uniquify=False):
                i, j = match[0], match[1]
                protected.add((min(i, j), max(i, j)))

        return protected

    @staticmethod
    def _fused_ring_atoms_to_remove(
        rdmol: Chem.Mol,
        torsion_region: Set[int],
    ) -> Set[int]:
        """Atoms on distant fused rings to remove.

        Fused-ring atoms not sharing a ring with any torsion-region
        atom are marked for removal, but only when the fused bond
        is C-C.
        """
        ring_info = rdmol.GetRingInfo()
        to_remove: Set[int] = set()

        # Identify atoms in multiple rings (fused atoms) and in any ring.
        atoms_in_ring: Set[int] = set()
        atoms_in_fused: Set[int] = set()
        for idx in range(rdmol.GetNumAtoms()):
            n_rings = RingInfo.NumAtomRings(ring_info, idx)
            if n_rings > 0:
                atoms_in_ring.add(idx)
            if n_rings > 1:
                atoms_in_fused.add(idx)

        # Only keep C-C bonds in fused ring positions as cutting candidates.
        fused_cc_bonds: List[Tuple[int, int]] = []
        for bond in rdmol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if i in atoms_in_fused and j in atoms_in_fused:
                if (
                    rdmol.GetAtomWithIdx(i).GetAtomicNum() == 6
                    and rdmol.GetAtomWithIdx(j).GetAtomicNum() == 6
                ):
                    fused_cc_bonds.append((i, j))

        # Collect all atoms in rings containing a fused C-C bond.
        all_atoms_in_fused_ring: Set[int] = set()
        for idx in atoms_in_ring:
            for bond_i, bond_j in fused_cc_bonds:
                for fused_idx in (bond_i, bond_j):
                    if RingInfo.AreAtomsInSameRing(ring_info, idx, fused_idx):
                        all_atoms_in_fused_ring.add(idx)

        # Remove fused ring atoms not in the same ring as torsion atoms.
        for idx in all_atoms_in_fused_ring:
            shares_ring = False
            for torsion_idx in torsion_region:
                if torsion_idx in all_atoms_in_fused_ring:
                    if RingInfo.AreAtomsInSameRing(ring_info, idx, torsion_idx):
                        shares_ring = True
                        break
            if not shares_ring:
                to_remove.add(idx)
                # Also remove attached hydrogens.
                for neig in rdmol.GetAtomWithIdx(idx).GetNeighbors():
                    if neig.GetAtomicNum() == 1:
                        to_remove.add(neig.GetIdx())

        return to_remove

    @staticmethod
    def _connected_components(
        adj: Dict[int, Set[int]],
        n_atoms: int,
    ) -> List[Set[int]]:
        """Return connected components as a list of atom-index sets."""
        visited: Set[int] = set()
        components: List[Set[int]] = []
        for start in range(n_atoms):
            if start in visited:
                continue
            # BFS
            comp: Set[int] = set()
            queue = [start]
            while queue:
                node = queue.pop()
                if node in visited:
                    continue
                visited.add(node)
                comp.add(node)
                for neig in adj.get(node, set()):
                    if neig not in visited:
                        queue.append(neig)
            if comp:
                components.append(comp)
        return components

    def __repr__(self) -> str:
        return (
            f"RuleBasedFragmenter("
            f"max_fragment_atoms={self._max_fragment_atoms})"
        )
