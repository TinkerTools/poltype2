"""
poltype.pipeline.stages.finalization – final output assembly stage.

Collects the key pipeline artifacts (optimised molecule, ESP results,
multipole frames, torsion scans) and assembles a summary of the
pipeline outputs.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger(__name__)


class FinalizationStage(Stage):
    """Assemble final output files (.xyz, .key).

    Gathers the artifacts produced by upstream stages and produces a
    ``final_summary`` artifact containing the collected outputs and
    output directory information.
    """

    def __init__(self) -> None:
        super().__init__(name="finalization")

    def execute(self, context: PipelineContext) -> StageResult:
        """Assemble final outputs from upstream artifacts.

        Returns
        -------
        StageResult
            ``COMPLETED`` with a ``final_summary`` artifact.
        """
        summary: Dict[str, Any] = {
            "work_dir": str(context.work_dir),
            "has_molecule": context.molecule is not None,
        }

        if context.molecule is not None:
            summary["molecule_name"] = context.molecule.name
            summary["num_atoms"] = context.molecule.num_atoms
            summary["formula"] = context.molecule.formula

        collected: List[str] = []

        opt_result = context.get_artifact("opt_result")
        if opt_result is not None:
            summary["opt_energy"] = opt_result.energy
            summary["opt_converged"] = opt_result.converged
            collected.append("opt_result")

        esp_result = context.get_artifact("esp_result")
        if esp_result is not None:
            summary["esp_grid_points"] = len(esp_result.esp_values)
            collected.append("esp_result")

        multipole_frames = context.get_artifact("multipole_frames")
        if multipole_frames is not None:
            summary["multipole_frames"] = multipole_frames
            collected.append("multipole_frames")

        torsion_scans = context.get_artifact("torsion_scans")
        if torsion_scans is not None:
            summary["torsion_bonds_scanned"] = len(torsion_scans)
            collected.append("torsion_scans")

        summary["collected_artifacts"] = collected

        logger.info(
            "Finalization complete: collected %d artifact(s): %s",
            len(collected),
            ", ".join(collected) if collected else "(none)",
        )

        return StageResult(
            status=StageStatus.COMPLETED,
            message=f"Finalized with {len(collected)} artifact(s)",
            artifacts={"final_summary": summary},
        )
