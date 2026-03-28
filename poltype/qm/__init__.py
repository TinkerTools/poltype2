"""
poltype.qm – QM backend abstraction and result types.

Re-exports the :class:`QMBackend` ABC, result dataclasses, and the
``select_backend`` factory so callers can write::

    from poltype.qm import QMBackend, select_backend
"""

from __future__ import annotations

from poltype.qm.backend import (
    ESPGridResult,
    OptimizationResult,
    QMBackend,
    TorsionScanResult,
    WBOResult,
)
from poltype.qm.factory import select_backend

__all__ = [
    "QMBackend",
    "OptimizationResult",
    "ESPGridResult",
    "TorsionScanResult",
    "WBOResult",
    "select_backend",
]
