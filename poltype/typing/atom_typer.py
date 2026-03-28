"""
poltype.typing.atom_typer – atom type assignment strategies.

:class:`AtomTyper` is the abstract interface.  Concrete implementations
match atoms to Tinker atom types.  :class:`SMARTSAtomTyper` reads an
ordered SMARTS→class mapping file and assigns types by last-match-wins.

Mapping file format (``amoeba21smartstoclass.txt``)
---------------------------------------------------
Each non-blank, non-comment line contains whitespace-separated tokens::

    SMARTS  classNumber  className  # optional comments

Example::

    [CX4H3]  35  C33  # C33, sp3 carbon with 3 Hydrogens
    [H][C]   16  HC   # HC, H on carbon

Lines starting with ``#`` are ignored.  Inline ``#`` starts a comment.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from rdkit import Chem

from poltype.typing.type_record import TypeRecord

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


class _SMARTSRule:
    """One line from a SMARTS mapping file."""

    __slots__ = ("smarts", "pattern", "atom_type", "atom_class", "description")

    def __init__(
        self,
        smarts: str,
        atom_type: int,
        atom_class: int,
        description: str,
    ) -> None:
        self.smarts = smarts
        self.pattern: Optional[Chem.Mol] = Chem.MolFromSmarts(smarts)
        self.atom_type = atom_type
        self.atom_class = atom_class
        self.description = description

    def __repr__(self) -> str:
        return (
            f"_SMARTSRule({self.smarts!r}, type={self.atom_type}, "
            f"class={self.atom_class})"
        )


# ---------------------------------------------------------------------------
# AtomTyper ABC
# ---------------------------------------------------------------------------


class AtomTyper(ABC):
    """Abstract base for atom type assignment strategies."""

    @abstractmethod
    def assign_types(self, rdmol: Chem.Mol) -> List[TypeRecord]:
        """Assign types to every atom in *rdmol*.

        Parameters
        ----------
        rdmol:
            RDKit molecule (should have explicit hydrogens).

        Returns
        -------
        list[TypeRecord]
            One record per atom, ordered by 0-based atom index.
        """


# ---------------------------------------------------------------------------
# SMARTSAtomTyper
# ---------------------------------------------------------------------------


class SMARTSAtomTyper(AtomTyper):
    """Assign atom types by matching SMARTS patterns.

    Rules are tested in order.  The **last** matching rule wins for each
    atom (later, more specific patterns override earlier, generic ones).
    Atoms that match no rule receive a fallback type of ``0`` with
    class ``0`` and description ``"Untyped"``.

    Parameters
    ----------
    rules:
        Ordered sequence of ``_SMARTSRule`` objects.  Normally populated
        via :meth:`from_file`.
    """

    def __init__(self, rules: Sequence[_SMARTSRule]) -> None:
        self._rules: List[_SMARTSRule] = list(rules)

    # -- alternate constructors --

    @classmethod
    def from_file(cls, path: str | Path) -> "SMARTSAtomTyper":
        """Load rules from the AMOEBA SMARTS-to-class mapping file.

        The expected format per line is::

            SMARTS  classNumber  className  # optional comment

        Parameters
        ----------
        path:
            Path to the mapping file.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"SMARTS mapping file not found: {path}")

        rules: List[_SMARTSRule] = []
        for lineno, raw_line in enumerate(path.read_text().splitlines(), 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            # Strip inline comments (after a # preceded by whitespace).
            comment_idx = line.find(" #")
            if comment_idx != -1:
                line = line[:comment_idx].strip()
            tokens = line.split()
            if len(tokens) < 2:
                continue
            smarts = tokens[0]
            try:
                class_num = int(tokens[1])
            except ValueError:
                logger.warning(
                    "Skipping line %d: class number is not an integer: %r",
                    lineno,
                    tokens[1],
                )
                continue
            class_name = tokens[2] if len(tokens) > 2 else ""
            rule = _SMARTSRule(smarts, class_num, class_num, class_name)
            if rule.pattern is None:
                logger.warning(
                    "Skipping invalid SMARTS %r at %s:%d",
                    smarts,
                    path,
                    lineno,
                )
                continue
            rules.append(rule)

        logger.debug("Loaded %d SMARTS rules from %s", len(rules), path)
        return cls(rules)

    @classmethod
    def from_lines(cls, lines: Sequence[str]) -> "SMARTSAtomTyper":
        """Build a typer from in-memory text lines (useful for testing).

        Accepts either the 3-column AMOEBA format
        (``SMARTS classNum className``) or the simpler 3-integer format
        (``SMARTS type class [description]``).
        """
        rules: List[_SMARTSRule] = []
        for raw_line in lines:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            # Strip inline comments.
            comment_idx = line.find(" #")
            if comment_idx != -1:
                line = line[:comment_idx].strip()
            tokens = line.split()
            if len(tokens) < 2:
                continue
            smarts = tokens[0]
            try:
                val1 = int(tokens[1])
            except ValueError:
                continue
            # If 3rd token is an int, treat as (type, class); else (class, name).
            atom_type = val1
            atom_class = val1
            description = ""
            if len(tokens) > 2:
                try:
                    atom_class = int(tokens[2])
                    description = " ".join(tokens[3:]) if len(tokens) > 3 else ""
                except ValueError:
                    # 3rd token is a class name string.
                    description = tokens[2]
            rule = _SMARTSRule(smarts, atom_type, atom_class, description)
            if rule.pattern is not None:
                rules.append(rule)
        return cls(rules)

    # -- properties --

    @property
    def num_rules(self) -> int:
        """Number of loaded SMARTS rules."""
        return len(self._rules)

    # -- core logic --

    def assign_types(self, rdmol: Chem.Mol) -> List[TypeRecord]:
        """Assign types by SMARTS scanning.

        Rules are applied in order.  Later (more specific) rules
        **override** earlier (generic) ones, matching the convention
        used by the legacy Poltype2 code.

        Unmatched atoms receive a fallback
        ``TypeRecord(atom_type=0, atom_class=0, description="Untyped")``.
        """
        n_atoms = rdmol.GetNumAtoms()
        assigned: Dict[int, TypeRecord] = {}

        for rule in self._rules:
            if rule.pattern is None:
                continue
            matches = rdmol.GetSubstructMatches(rule.pattern)
            for match_tuple in matches:
                for atom_idx in match_tuple:
                    # Last match wins (more specific rules come later).
                    assigned[atom_idx] = TypeRecord(
                        atom_index=atom_idx,
                        atom_type=rule.atom_type,
                        atom_class=rule.atom_class,
                        description=rule.description,
                        smarts=rule.smarts,
                    )

        # Fill unmatched atoms with a fallback record.
        records: List[TypeRecord] = []
        for i in range(n_atoms):
            if i in assigned:
                records.append(assigned[i])
            else:
                records.append(
                    TypeRecord(
                        atom_index=i,
                        atom_type=0,
                        atom_class=0,
                        description="Untyped",
                    )
                )
        return records

    def __repr__(self) -> str:
        return f"SMARTSAtomTyper(rules={self.num_rules})"
