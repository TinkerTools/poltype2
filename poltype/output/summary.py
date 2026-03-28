"""
poltype.output.summary – pipeline run summary generation.

Produces a structured :class:`RunSummary` from the pipeline context and
writes a human-readable text report via :class:`SummaryWriter`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List

if TYPE_CHECKING:
    from poltype.pipeline.context import PipelineContext

from poltype.output.writer import OutputWriter

logger = logging.getLogger(__name__)


@dataclass
class RunSummary:
    """Structured summary of a pipeline run.

    Attributes
    ----------
    molecule_name:
        Name of the input molecule (or ``"N/A"``).
    formula:
        Molecular formula (or ``"N/A"``).
    num_atoms:
        Number of atoms.
    backend_name:
        Name of the QM backend used (or ``"none"``).
    work_dir:
        Working directory path.
    stage_results:
        Mapping of stage name → (status, message, elapsed_seconds).
    output_files:
        Paths to written output files.
    timestamp:
        ISO-8601 timestamp of summary generation.
    """

    molecule_name: str = "N/A"
    formula: str = "N/A"
    num_atoms: int = 0
    backend_name: str = "none"
    work_dir: str = ""
    stage_results: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    output_files: List[str] = field(default_factory=list)
    timestamp: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(),
    )

    @classmethod
    def from_context(
        cls,
        context: "PipelineContext",
        output_files: List[Path] | None = None,
    ) -> RunSummary:
        """Build a summary from a completed pipeline context.

        Parameters
        ----------
        context:
            Pipeline context after execution.
        output_files:
            Paths to any files written by output writers.
        """
        mol = context.molecule
        stage_results: Dict[str, Dict[str, Any]] = {}
        for name, result in context.stage_results.items():
            stage_results[name] = {
                "status": result.status.value,
                "message": result.message,
                "elapsed_seconds": round(result.elapsed_seconds, 3),
            }

        return cls(
            molecule_name=mol.name if mol else "N/A",
            formula=mol.formula if mol else "N/A",
            num_atoms=mol.num_atoms if mol else 0,
            backend_name=context.backend.name if context.backend else "none",
            work_dir=str(context.work_dir),
            stage_results=stage_results,
            output_files=[str(p) for p in (output_files or [])],
        )


class SummaryWriter(OutputWriter):
    """Write a human-readable pipeline summary to a text file.

    Parameters
    ----------
    output_dir:
        Directory for the output file.  Falls back to ``context.work_dir``.
    filename:
        Name of the summary file.  Defaults to ``"poltype_summary.txt"``.
    """

    def __init__(
        self,
        output_dir: Path | None = None,
        filename: str = "poltype_summary.txt",
    ) -> None:
        super().__init__(output_dir=output_dir)
        self._filename = filename
        self._extra_output_files: List[Path] = []

    @property
    def name(self) -> str:
        return "summary_writer"

    def set_output_files(self, files: List[Path]) -> None:
        """Record output file paths written by other writers."""
        self._extra_output_files = list(files)

    def write(self, context: "PipelineContext") -> Path:
        """Write the summary text file.

        Returns
        -------
        Path
            Absolute path to the written summary file.
        """
        out_dir = self._resolve_dir(context)
        out_path = out_dir / self._filename

        summary = RunSummary.from_context(
            context, output_files=self._extra_output_files,
        )
        text = self._format(summary)
        out_path.write_text(text)
        logger.info("Wrote pipeline summary: %s", out_path)
        return out_path

    def _format(self, summary: RunSummary) -> str:
        """Format the summary as a human-readable string."""
        sep = "=" * 60
        lines = [
            sep,
            "  Poltype2 Pipeline Run Summary",
            sep,
            "",
            f"  Timestamp : {summary.timestamp}",
            f"  Molecule  : {summary.molecule_name}",
            f"  Formula   : {summary.formula}",
            f"  Atoms     : {summary.num_atoms}",
            f"  Backend   : {summary.backend_name}",
            f"  Work dir  : {summary.work_dir}",
            "",
            "-" * 60,
            "  Stage Results",
            "-" * 60,
        ]

        for stage_name, info in summary.stage_results.items():
            status = info["status"].upper()
            elapsed = info["elapsed_seconds"]
            message = info["message"]
            lines.append(
                f"  {stage_name:<30s}  {status:<10s}  "
                f"{elapsed:>7.3f}s  {message}"
            )

        lines.append("")

        if summary.output_files:
            lines.append("-" * 60)
            lines.append("  Output Files")
            lines.append("-" * 60)
            for f in summary.output_files:
                lines.append(f"  {f}")
            lines.append("")

        lines.append(sep)
        lines.append("")
        return "\n".join(lines)
