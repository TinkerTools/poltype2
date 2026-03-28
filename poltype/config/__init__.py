"""
poltype.config – configuration schema and loading utilities.

Re-exports the typed dataclasses and the ``load_config`` helper so
callers can write::

    from poltype.config import PoltypeConfig, load_config
"""

from __future__ import annotations

from poltype.config.loader import load_config
from poltype.config.schema import PoltypeConfig, QMConfig, ResourceConfig

__all__ = ["PoltypeConfig", "QMConfig", "ResourceConfig", "load_config"]
