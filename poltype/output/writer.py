"""
poltype.output.writer – abstract base class for output writers.

Every output writer inherits from :class:`OutputWriter` and implements
the ``write`` method.  Writers receive the pipeline context (including
all artifacts) and produce a single output file.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from poltype.pipeline.context import PipelineContext


class OutputWriter(ABC):
    """Abstract base class for pipeline output writers.

    Subclasses **must** implement :meth:`write` and the :attr:`name`
    property.

    Parameters
    ----------
    output_dir:
        Directory where output files will be written.  If ``None``,
        each writer should use the pipeline's ``work_dir``.
    """

    def __init__(self, output_dir: Path | None = None) -> None:
        self._output_dir = output_dir

    @property
    def output_dir(self) -> Path | None:
        """Target directory for written files."""
        return self._output_dir

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable name of this writer (e.g. ``"key_writer"``)."""

    @abstractmethod
    def write(self, context: "PipelineContext") -> Path:
        """Write an output file from the pipeline context.

        Parameters
        ----------
        context:
            Pipeline context containing config, molecule, backend, and
            inter-stage artifacts.

        Returns
        -------
        Path
            Absolute path to the written file.
        """

    def _resolve_dir(self, context: "PipelineContext") -> Path:
        """Return the output directory, falling back to ``context.work_dir``."""
        d = self._output_dir if self._output_dir is not None else context.work_dir
        d.mkdir(parents=True, exist_ok=True)
        return d

    def __repr__(self) -> str:
        return f"{type(self).__name__}(output_dir={self._output_dir})"
