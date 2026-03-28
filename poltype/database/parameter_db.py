"""
poltype.database.parameter_db – parameter database loading & lookup.

:class:`ParameterDatabase` is the abstract interface.
:class:`SmallMoleculeDB` loads the AMOEBA SMARTS→type mapping shipped
with Poltype2 and reports which atoms have existing parameters.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Optional

from rdkit import Chem

from poltype.database.match_result import DatabaseMatchResult
from poltype.typing.atom_typer import SMARTSAtomTyper
from poltype.typing.type_record import TypeRecord

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Abstract base
# ---------------------------------------------------------------------------


class ParameterDatabase(ABC):
    """Abstract interface for parameter databases."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable database identifier."""

    @abstractmethod
    def lookup(
        self,
        rdmol: Chem.Mol,
        type_records: Optional[List[TypeRecord]] = None,
    ) -> DatabaseMatchResult:
        """Look up which atoms/terms already have parameters.

        Parameters
        ----------
        rdmol:
            RDKit molecule.
        type_records:
            Optional pre-computed type records.  If provided the database
            can use them to accelerate or refine the search.

        Returns
        -------
        DatabaseMatchResult
        """


# ---------------------------------------------------------------------------
# Small-molecule database backed by SMARTS mapping file
# ---------------------------------------------------------------------------


class SmallMoleculeDB(ParameterDatabase):
    """AMOEBA small-molecule parameter database.

    Wraps a SMARTS→type mapping file and reports which atoms have a
    matching type in the database (i.e. have pre-existing parameters).

    Parameters
    ----------
    smarts_map_path:
        Path to the SMARTS→type mapping file
        (e.g. ``amoeba21smartstoclass.txt``).
    """

    def __init__(self, smarts_map_path: str | Path) -> None:
        self._path = Path(smarts_map_path)
        self._typer = SMARTSAtomTyper.from_file(self._path)

    @property
    def name(self) -> str:
        return f"SmallMoleculeDB({self._path.name})"

    @property
    def num_rules(self) -> int:
        """Number of rules loaded from the mapping file."""
        return self._typer.num_rules

    def lookup(
        self,
        rdmol: Chem.Mol,
        type_records: Optional[List[TypeRecord]] = None,
    ) -> DatabaseMatchResult:
        """Determine which atoms have a known type in the database.

        Atoms whose assigned type is non-zero (i.e. matched a SMARTS
        rule) are reported as *matched*; all others as *unmatched*.
        """
        records = type_records or self._typer.assign_types(rdmol)
        matched: List[int] = []
        unmatched: List[int] = []
        for rec in records:
            if rec.atom_type != 0:
                matched.append(rec.atom_index)
            else:
                unmatched.append(rec.atom_index)

        result = DatabaseMatchResult(
            matched_atom_indices=matched,
            unmatched_atom_indices=unmatched,
            source=self.name,
            metadata={"num_rules": self._typer.num_rules},
        )
        logger.debug(
            "Database lookup: %d/%d atoms matched (%s)",
            len(matched),
            len(matched) + len(unmatched),
            self.name,
        )
        return result

    def __repr__(self) -> str:
        return f"SmallMoleculeDB(path={self._path}, rules={self.num_rules})"
