"""
poltype.pipeline.events – pipeline event / hook system.

The event system lets external code observe pipeline execution without
coupling to the runner internals.  An :class:`EventBus` holds an ordered
list of :class:`EventHook` callables and broadcasts :class:`EventData`
objects to them.

Usage::

    from poltype.pipeline.events import EventBus, PipelineEvent

    bus = EventBus()
    bus.register(lambda data: print(data.event, data.stage_name))
"""

from __future__ import annotations

import enum
import logging
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

logger = logging.getLogger("poltype.pipeline.events")


# ---------------------------------------------------------------------------
# Event types
# ---------------------------------------------------------------------------


class PipelineEvent(enum.Enum):
    """Events emitted during pipeline execution."""

    PIPELINE_STARTED = "pipeline_started"
    PIPELINE_COMPLETED = "pipeline_completed"
    STAGE_STARTED = "stage_started"
    STAGE_COMPLETED = "stage_completed"
    STAGE_SKIPPED = "stage_skipped"
    STAGE_FAILED = "stage_failed"


# ---------------------------------------------------------------------------
# Event data
# ---------------------------------------------------------------------------


@dataclass
class EventData:
    """Payload delivered to every registered hook.

    Attributes
    ----------
    event:
        The event that occurred.
    stage_name:
        Name of the stage that triggered the event (empty string for
        pipeline-level events).
    message:
        Human-readable description.
    timestamp:
        ``time.monotonic()`` when the event was created.
    extra:
        Arbitrary key→value metadata (e.g. elapsed seconds, result).
    """

    event: PipelineEvent
    stage_name: str = ""
    message: str = ""
    timestamp: float = 0.0
    extra: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.timestamp == 0.0:
            self.timestamp = time.monotonic()


# ---------------------------------------------------------------------------
# Hook type
# ---------------------------------------------------------------------------

#: A hook is any callable that accepts an :class:`EventData` instance.
EventHook = Callable[[EventData], None]


# ---------------------------------------------------------------------------
# Event bus
# ---------------------------------------------------------------------------


class EventBus:
    """Registry of :data:`EventHook` callables.

    Hooks are invoked in the order they were registered.  A hook that
    raises an exception is logged and skipped — it does **not** abort the
    pipeline.

    Parameters
    ----------
    hooks:
        Optional initial list of hooks.
    """

    def __init__(self, hooks: Optional[List[EventHook]] = None) -> None:
        self._hooks: List[EventHook] = list(hooks) if hooks else []

    def register(self, hook: EventHook) -> None:
        """Register a new hook."""
        self._hooks.append(hook)

    def emit(self, data: EventData) -> None:
        """Broadcast *data* to all registered hooks.

        Exceptions raised by individual hooks are caught and logged so
        that a misbehaving hook cannot crash the pipeline.
        """
        for hook in self._hooks:
            try:
                hook(data)
            except Exception:
                logger.exception(
                    "Hook %r raised an exception for event %s",
                    hook,
                    data.event,
                )

    @property
    def hooks(self) -> List[EventHook]:
        """Read-only copy of the registered hooks."""
        return list(self._hooks)

    def __repr__(self) -> str:
        return f"EventBus(hooks={len(self._hooks)})"
