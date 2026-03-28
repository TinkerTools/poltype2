"""
poltype.pipeline.stages.atom_typing – atom type assignment stage.

Assigns Tinker atom types and classes to every atom in the current
molecule using a SMARTS-based matching strategy.  The resulting
:class:`TypeRecord` list is stored as the ``"type_records"`` artifact
for consumption by downstream stages (database matching, output
writers).

This stage is skipped when ``config.database_match_only`` is ``True``
*and* pre-computed type records are already present in the context, or
when the pipeline is running in dry-run mode and no molecule is
available.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, List, Optional

from poltype.pipeline.stage import Stage, StageResult, StageStatus
from poltype.typing.atom_typer import AtomTyper, SMARTSAtomTyper
from poltype.typing.type_record import TypeRecord

if TYPE_CHECKING:
    from poltype.pipeline.context import PipelineContext

logger = logging.getLogger(__name__)


class AtomTypingStage(Stage):
    """Assign atom types via SMARTS matching.

    Parameters
    ----------
    typer:
        An :class:`AtomTyper` instance.  If ``None``, one is built
        automatically from ``config.small_molecule_smarts_to_tinker_class``
        during execution.
    """

    def __init__(self, typer: Optional[AtomTyper] = None) -> None:
        super().__init__(name="atom_typing")
        self._typer = typer

    def execute(self, context: "PipelineContext") -> StageResult:
        mol = context.molecule
        if mol is None:
            return StageResult(
                status=StageStatus.FAILED,
                message="No molecule in context",
            )

        # Build typer lazily from config if not injected.
        typer = self._typer
        if typer is None:
            smarts_path = context.config.small_molecule_smarts_to_tinker_class
            try:
                typer = SMARTSAtomTyper.from_file(smarts_path)
            except FileNotFoundError:
                logger.warning(
                    "SMARTS mapping file not found: %s; "
                    "falling back to index-based types",
                    smarts_path,
                )
                # Produce fallback records (type = tinker index).
                records = _fallback_records(mol.rdmol)
                context.set_artifact("type_records", records)
                return StageResult(
                    status=StageStatus.COMPLETED,
                    message=f"Fallback typing for {mol.num_atoms} atoms",
                    artifacts={"type_records": records},
                )

        records: List[TypeRecord] = typer.assign_types(mol.rdmol)
        context.set_artifact("type_records", records)

        n_typed = sum(1 for r in records if r.atom_type != 0)
        logger.info(
            "Atom typing complete: %d/%d atoms matched",
            n_typed,
            mol.num_atoms,
        )
        return StageResult(
            status=StageStatus.COMPLETED,
            message=f"Typed {n_typed}/{mol.num_atoms} atoms",
            artifacts={"type_records": records},
        )

    def should_skip(self, context: "PipelineContext") -> bool:
        """Skip in dry-run mode."""
        return getattr(context.config, "dry_run", False)


def _fallback_records(rdmol) -> List[TypeRecord]:
    """Generate index-based fallback type records."""
    records = []
    for i in range(rdmol.GetNumAtoms()):
        atom = rdmol.GetAtomWithIdx(i)
        tinker_idx = i + 1
        records.append(
            TypeRecord(
                atom_index=i,
                atom_type=tinker_idx,
                atom_class=tinker_idx,
                description=f"{atom.GetSymbol()}{tinker_idx}",
            )
        )
    return records
