"""
poltype.pipeline.context – mutable per-run pipeline state.

A :class:`PipelineContext` is created once at the beginning of a
pipeline run and flows through every stage.  It carries:

* the immutable :class:`PoltypeConfig`,
* the current :class:`Molecule` (which may be replaced after geometry
  optimisation),
* the selected :class:`QMBackend`,
* an ``artifacts`` dict for inter-stage data exchange, and
* a ``stage_results`` dict that records the outcome of each stage.

This module has **no dependency on PoltypeModules** and can be
unit-tested in isolation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from poltype.config.schema import PoltypeConfig
from poltype.molecule.molecule import Molecule
from poltype.pipeline.stage import StageResult
from poltype.qm.backend import QMBackend


class PipelineContext:
    """Mutable per-run state that flows through every pipeline stage.

    Parameters
    ----------
    config:
        Immutable run configuration.
    molecule:
        Initial molecule (may be ``None`` if loaded during the
        input-preparation stage).
    backend:
        QM backend to use for calculations (may be ``None`` if
        selected during the pipeline).
    work_dir:
        Working directory for file I/O.  Defaults to the current
        working directory.
    """

    def __init__(
        self,
        config: PoltypeConfig,
        molecule: Optional[Molecule] = None,
        backend: Optional[QMBackend] = None,
        work_dir: Optional[Path] = None,
    ) -> None:
        self._config = config
        self._molecule = molecule
        self._backend = backend
        self._work_dir = Path(work_dir).resolve() if work_dir else Path.cwd()
        self.artifacts: Dict[str, Any] = {}
        self.stage_results: Dict[str, StageResult] = {}

    # -- read-only property for config --

    @property
    def config(self) -> PoltypeConfig:
        """Immutable run configuration."""
        return self._config

    # -- molecule (reassignable) --

    @property
    def molecule(self) -> Optional[Molecule]:
        """Current molecule (may be replaced after geometry optimisation)."""
        return self._molecule

    @molecule.setter
    def molecule(self, value: Optional[Molecule]) -> None:
        self._molecule = value

    # -- backend (reassignable) --

    @property
    def backend(self) -> Optional[QMBackend]:
        """QM backend for the run."""
        return self._backend

    @backend.setter
    def backend(self, value: Optional[QMBackend]) -> None:
        self._backend = value

    # -- work_dir --

    @property
    def work_dir(self) -> Path:
        """Working directory for file I/O."""
        return self._work_dir

    # -- artifact helpers --

    def get_artifact(self, key: str, default: Any = None) -> Any:
        """Retrieve an artifact by key, returning *default* if absent."""
        return self.artifacts.get(key, default)

    def set_artifact(self, key: str, value: Any) -> None:
        """Store an artifact for downstream stages."""
        self.artifacts[key] = value

    # -- molecule replacement --

    def update_molecule(self, new_mol: Molecule) -> None:
        """Replace the current molecule (e.g. after geometry optimisation).

        Parameters
        ----------
        new_mol:
            The new molecule to use for subsequent stages.
        """
        self._molecule = new_mol

    def __repr__(self) -> str:
        mol_name = self._molecule.name if self._molecule else "None"
        return (
            f"PipelineContext(config=PoltypeConfig, "
            f"molecule={mol_name!r}, "
            f"work_dir={str(self._work_dir)!r})"
        )
