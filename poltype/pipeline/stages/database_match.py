"""
poltype.pipeline.stages.database_match – parameter database lookup stage.

Looks up existing parameters for the molecule in the small-molecule
database.  The :class:`DatabaseMatchResult` is stored as the
``"db_match_result"`` artifact so downstream stages can decide whether
fitting is needed or database values can be re-used.

This stage is skipped in dry-run mode.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional

from poltype.database.match_result import DatabaseMatchResult
from poltype.database.parameter_db import ParameterDatabase, SmallMoleculeDB
from poltype.pipeline.stage import Stage, StageResult, StageStatus

if TYPE_CHECKING:
    from poltype.pipeline.context import PipelineContext

logger = logging.getLogger(__name__)


class DatabaseMatchStage(Stage):
    """Look up existing parameters in the AMOEBA database.

    Parameters
    ----------
    database:
        A :class:`ParameterDatabase` instance.  If ``None``, a
        :class:`SmallMoleculeDB` is built from the config path
        ``small_molecule_smarts_to_tinker_class`` during execution.
    """

    def __init__(self, database: Optional[ParameterDatabase] = None) -> None:
        super().__init__(name="database_match")
        self._database = database

    def execute(self, context: "PipelineContext") -> StageResult:
        mol = context.molecule
        if mol is None:
            return StageResult(
                status=StageStatus.FAILED,
                message="No molecule in context",
            )

        # Build database lazily from config if not injected.
        db = self._database
        if db is None:
            smarts_path = context.config.small_molecule_smarts_to_tinker_class
            try:
                db = SmallMoleculeDB(smarts_path)
            except FileNotFoundError:
                logger.warning(
                    "Database file not found: %s; skipping lookup",
                    smarts_path,
                )
                result = DatabaseMatchResult(
                    matched_atom_indices=[],
                    unmatched_atom_indices=list(range(mol.num_atoms)),
                    source="none",
                )
                context.set_artifact("db_match_result", result)
                return StageResult(
                    status=StageStatus.COMPLETED,
                    message="No database available; all atoms unmatched",
                    artifacts={"db_match_result": result},
                )

        # Use pre-computed type records if available.
        type_records = context.get_artifact("type_records")
        result: DatabaseMatchResult = db.lookup(mol.rdmol, type_records)
        context.set_artifact("db_match_result", result)

        logger.info(
            "Database match: %.0f%% coverage (%s)",
            result.coverage * 100,
            result.source,
        )
        return StageResult(
            status=StageStatus.COMPLETED,
            message=(
                f"{len(result.matched_atom_indices)}/{mol.num_atoms} atoms "
                f"matched ({result.source})"
            ),
            artifacts={"db_match_result": result},
        )

    def should_skip(self, context: "PipelineContext") -> bool:
        """Skip in dry-run mode."""
        return getattr(context.config, "dry_run", False)
