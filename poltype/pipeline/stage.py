"""
poltype.pipeline.stage – abstract stage interface and result types.

Every pipeline stage inherits from :class:`Stage` and implements the
``execute`` method.  The :class:`StageResult` dataclass carries the
outcome of a single stage execution, including an ``artifacts`` dict
that downstream stages can consume.
"""

from __future__ import annotations

import enum
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict

if TYPE_CHECKING:
    from poltype.pipeline.context import PipelineContext


# ---------------------------------------------------------------------------
# Stage status
# ---------------------------------------------------------------------------


class StageStatus(enum.Enum):
    """Lifecycle status of a pipeline stage."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    SKIPPED = "skipped"
    FAILED = "failed"


# ---------------------------------------------------------------------------
# Stage result
# ---------------------------------------------------------------------------


@dataclass
class StageResult:
    """Outcome of executing a single pipeline stage.

    Attributes
    ----------
    status:
        Final status after execution.
    message:
        Human-readable note (e.g. error description, skip reason).
    elapsed_seconds:
        Wall-clock time the stage took to execute.
    artifacts:
        Arbitrary key→value map of outputs the stage wants to expose
        to downstream stages via :pyattr:`PipelineContext.artifacts`.
    """

    status: StageStatus
    message: str = ""
    elapsed_seconds: float = 0.0
    artifacts: Dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Stage ABC
# ---------------------------------------------------------------------------


class Stage(ABC):
    """Abstract base class for a composable pipeline stage.

    Subclasses **must** implement :meth:`execute`.  They **may** override
    :meth:`should_skip` to declare config-driven skip conditions.

    Parameters
    ----------
    name:
        Unique name for this stage (used as a key in
        ``PipelineContext.stage_results``).
    """

    def __init__(self, name: str) -> None:
        self._name = name

    @property
    def name(self) -> str:
        """Unique stage identifier."""
        return self._name

    @abstractmethod
    def execute(self, context: PipelineContext) -> StageResult:
        """Run the stage logic.

        Parameters
        ----------
        context:
            Mutable pipeline context carrying config, molecule,
            backend, and inter-stage artifacts.

        Returns
        -------
        StageResult
            Outcome including status and any artifacts to publish.
        """

    def should_skip(self, context: PipelineContext) -> bool:
        """Return ``True`` if this stage should be skipped.

        The default implementation always returns ``False``.  Subclasses
        can override to inspect ``context.config`` flags.
        """
        return False

    def __repr__(self) -> str:
        return f"Stage(name={self._name!r})"
