"""
poltype.pipeline.stages.finalization – final output assembly stage.

Collects the key pipeline artifacts (optimised molecule, ESP results,
multipole frames, torsion scans) and assembles a summary of the
pipeline outputs.  When output writers are configured, writes the
actual output files (e.g. ``.xyz``, ``.key``, summary report).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from poltype.output.writer import OutputWriter

from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger(__name__)


class FinalizationStage(Stage):
    """Assemble final output files (.xyz, .key) and a run summary.

    Gathers the artifacts produced by upstream stages, runs any
    configured :class:`OutputWriter` instances to produce actual files,
    and produces a ``final_summary`` artifact containing the collected
    outputs and output directory information.

    Parameters
    ----------
    writers:
        Optional list of output writers.  If ``None``, the stage uses
        the default set: :class:`XYZFileWriter`, :class:`KeyFileWriter`,
        and :class:`SummaryWriter`.
    """

    def __init__(
        self,
        writers: Optional[List[OutputWriter]] = None,
    ) -> None:
        super().__init__(name="finalization")
        self._writers: List[OutputWriter] = (
            list(writers) if writers is not None else []
        )

    @property
    def writers(self) -> List[OutputWriter]:
        """Configured output writers."""
        return list(self._writers)

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

        # --- run output writers ---
        written_files: List[Path] = []
        if self._writers and context.molecule is not None:
            for writer in self._writers:
                try:
                    path = writer.write(context)
                    written_files.append(path)
                    logger.info(
                        "Writer '%s' produced: %s", writer.name, path,
                    )
                except Exception as exc:
                    logger.warning(
                        "Writer '%s' failed: %s", writer.name, exc,
                    )

        summary["output_files"] = [str(p) for p in written_files]

        logger.info(
            "Finalization complete: collected %d artifact(s): %s, "
            "wrote %d file(s)",
            len(collected),
            ", ".join(collected) if collected else "(none)",
            len(written_files),
        )

        return StageResult(
            status=StageStatus.COMPLETED,
            message=(
                f"Finalized with {len(collected)} artifact(s), "
                f"{len(written_files)} file(s)"
            ),
            artifacts={
                "final_summary": summary,
                "output_files": written_files,
            },
        )
