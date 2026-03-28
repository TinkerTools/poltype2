"""
poltype.pipeline.checkpoint – save / restore pipeline progress.

A :class:`CheckpointManager` writes a small JSON manifest after every
successful stage so the pipeline can be **resumed** from the last
completed stage instead of re-running everything from scratch.

The checkpoint is a single file (``poltype_checkpoint.json``) written to
the working directory.  It records:

* which stages have completed (with their status and message),
* a monotonic timestamp, and
* the pipeline run ID (a random hex token).

Heavy artefacts (coordinates, ESP grids, …) are **not** serialised —
only the stage names and statuses.  On resume the runner simply skips
stages that are already recorded as ``COMPLETED`` in the checkpoint.

This module has **no dependency on PoltypeModules** and can be
unit-tested in isolation.
"""

from __future__ import annotations

import json
import logging
import secrets
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from poltype.pipeline.stage import StageResult, StageStatus

logger = logging.getLogger("poltype.pipeline.checkpoint")

#: Default checkpoint filename.
CHECKPOINT_FILENAME = "poltype_checkpoint.json"


# ---------------------------------------------------------------------------
# Checkpoint data
# ---------------------------------------------------------------------------


@dataclass
class CheckpointData:
    """Serialisable snapshot of pipeline progress.

    Attributes
    ----------
    run_id:
        Random hex token identifying this run.
    completed_stages:
        Ordered list of stage names that have completed successfully.
    stage_statuses:
        Mapping of stage name → status string for all recorded stages.
    stage_messages:
        Mapping of stage name → human-readable result message.
    timestamp:
        ``time.monotonic()`` when the checkpoint was written.
    """

    run_id: str = ""
    completed_stages: List[str] = field(default_factory=list)
    stage_statuses: Dict[str, str] = field(default_factory=dict)
    stage_messages: Dict[str, str] = field(default_factory=dict)
    timestamp: float = 0.0


# ---------------------------------------------------------------------------
# Checkpoint manager
# ---------------------------------------------------------------------------


class CheckpointManager:
    """Persists and restores pipeline progress.

    Parameters
    ----------
    checkpoint_dir:
        Directory where the checkpoint file is stored.  Defaults to the
        current working directory.
    filename:
        Name of the checkpoint file.  Defaults to
        ``poltype_checkpoint.json``.
    """

    def __init__(
        self,
        checkpoint_dir: Optional[Path] = None,
        filename: str = CHECKPOINT_FILENAME,
    ) -> None:
        self._dir = Path(checkpoint_dir).resolve() if checkpoint_dir else Path.cwd()
        self._filename = filename
        self._run_id: str = secrets.token_hex(8)
        self._data: CheckpointData = CheckpointData(run_id=self._run_id)

    @property
    def checkpoint_path(self) -> Path:
        """Full path to the checkpoint file."""
        return self._dir / self._filename

    @property
    def data(self) -> CheckpointData:
        """Current in-memory checkpoint data."""
        return self._data

    # -- persistence ---------------------------------------------------------

    def save(self, stage_name: str, result: StageResult) -> None:
        """Record a stage result and persist the checkpoint.

        Only ``COMPLETED`` stages are added to ``completed_stages``;
        all statuses are recorded in ``stage_statuses``.

        Parameters
        ----------
        stage_name:
            Name of the stage that just finished.
        result:
            The :class:`StageResult` returned by the stage.
        """
        self._data.stage_statuses[stage_name] = result.status.value
        self._data.stage_messages[stage_name] = result.message

        if result.status is StageStatus.COMPLETED:
            if stage_name not in self._data.completed_stages:
                self._data.completed_stages.append(stage_name)

        self._data.timestamp = time.monotonic()
        self._write()

    def load(self) -> CheckpointData:
        """Load a checkpoint from disk.

        Returns
        -------
        CheckpointData
            Restored checkpoint.  The internal ``_data`` attribute is
            also updated.

        Raises
        ------
        FileNotFoundError
            If the checkpoint file does not exist.
        """
        path = self.checkpoint_path
        raw = json.loads(path.read_text(encoding="utf-8"))
        self._data = CheckpointData(
            run_id=raw.get("run_id", ""),
            completed_stages=raw.get("completed_stages", []),
            stage_statuses=raw.get("stage_statuses", {}),
            stage_messages=raw.get("stage_messages", {}),
            timestamp=raw.get("timestamp", 0.0),
        )
        self._run_id = self._data.run_id
        logger.info(
            "Loaded checkpoint (run %s) with %d completed stages",
            self._run_id,
            len(self._data.completed_stages),
        )
        return self._data

    def should_skip_stage(self, stage_name: str) -> bool:
        """Return ``True`` if *stage_name* was already completed.

        This is used by the runner to skip stages that completed in a
        previous invocation.
        """
        return stage_name in self._data.completed_stages

    def clear(self) -> None:
        """Delete the checkpoint file (if it exists) and reset state."""
        path = self.checkpoint_path
        if path.exists():
            path.unlink()
            logger.info("Checkpoint file removed: %s", path)
        self._data = CheckpointData(run_id=self._run_id)

    # -- internal ------------------------------------------------------------

    def _write(self) -> None:
        """Serialise the current data to JSON."""
        payload = {
            "run_id": self._data.run_id,
            "completed_stages": self._data.completed_stages,
            "stage_statuses": self._data.stage_statuses,
            "stage_messages": self._data.stage_messages,
            "timestamp": self._data.timestamp,
        }
        self._dir.mkdir(parents=True, exist_ok=True)
        self.checkpoint_path.write_text(
            json.dumps(payload, indent=2) + "\n",
            encoding="utf-8",
        )
        logger.debug("Checkpoint saved to %s", self.checkpoint_path)

    def __repr__(self) -> str:
        n = len(self._data.completed_stages)
        return f"CheckpointManager(run={self._run_id!r}, completed={n})"
