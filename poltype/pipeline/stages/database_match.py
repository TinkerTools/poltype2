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
from poltype.output.params import MultipoleParam, PolarizeParam
from poltype.pipeline.parameter_utils import atom_type_for_index, torsion_indices_for_bond
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
        result.matched_torsions = self._matched_torsions(context, result)
        context.set_artifact("db_match_result", result)

        artifacts = {"db_match_result": result}
        multipoles = self._materialize_multipoles(context, type_records)
        if multipoles:
            artifacts["multipoles"] = multipoles

        polarization_params = self._materialize_polarization(context, type_records)
        if polarization_params:
            artifacts["polarization_params"] = polarization_params

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
            artifacts=artifacts,
        )

    def should_skip(self, context: "PipelineContext") -> bool:
        """Skip in dry-run mode."""
        return getattr(context.config, "dry_run", False)

    @staticmethod
    def _materialize_multipoles(
        context: "PipelineContext",
        type_records,
    ) -> list[MultipoleParam]:
        raw_records = context.get_artifact("fitted_multipoles_by_atom_index") or []
        if not raw_records:
            return []

        type_map = _type_map(type_records)
        multipoles_by_type: dict[int, MultipoleParam] = {}
        for record in raw_records:
            atom_type = atom_type_for_index(record["atom_index"], type_map)
            frame = [
                atom_type_for_index(frame_idx, type_map)
                for frame_idx in record["frame_atom_indices"]
            ]
            multipoles_by_type.setdefault(
                atom_type,
                MultipoleParam(
                    atom_type=atom_type,
                    frame=frame,
                    charge=record["charge"],
                    dipole=record["dipole"],
                    quadrupole=record["quadrupole"],
                ),
            )
        return sorted(multipoles_by_type.values(), key=lambda mp: mp.atom_type)

    @staticmethod
    def _materialize_polarization(
        context: "PipelineContext",
        type_records,
    ) -> list[PolarizeParam]:
        raw_records = context.get_artifact("polarization_params_by_atom_index") or []
        if not raw_records:
            return []

        type_map = _type_map(type_records)
        polarization_by_type: dict[int, PolarizeParam] = {}
        for record in raw_records:
            atom_type = atom_type_for_index(record["atom_index"], type_map)
            neighbor_types = sorted(
                {
                    atom_type_for_index(neighbor_idx, type_map)
                    for neighbor_idx in record["neighbor_atom_indices"]
                }
            )
            existing = polarization_by_type.get(atom_type)
            if existing is None:
                polarization_by_type[atom_type] = PolarizeParam(
                    atom_type=atom_type,
                    alpha=record["alpha"],
                    thole=record["thole"],
                    neighbors=neighbor_types,
                )
            else:
                existing.neighbors = sorted(set(existing.neighbors) | set(neighbor_types))
        return sorted(polarization_by_type.values(), key=lambda pp: pp.atom_type)

    @staticmethod
    def _matched_torsions(
        context: "PipelineContext",
        result: DatabaseMatchResult,
    ) -> set[tuple[int, int, int, int]]:
        """Infer which rotatable bonds can be re-used from database typing."""
        mol = context.molecule
        if mol is None:
            return set()

        matched_atoms = set(result.matched_atom_indices)
        matched_torsions: set[tuple[int, int, int, int]] = set(result.matched_torsions)
        for bond in mol.rotatable_bonds:
            torsion = torsion_indices_for_bond(mol.rdmol, bond)
            if torsion is None:
                continue
            if all(idx in matched_atoms for idx in torsion):
                matched_torsions.add(torsion)
        return matched_torsions


def _type_map(type_records) -> dict[int, int]:
    if not type_records:
        return {}
    return {record.atom_index: record.atom_type for record in type_records}
