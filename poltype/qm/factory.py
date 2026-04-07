"""
poltype.qm.factory – select the appropriate QM backend from config.

Usage::

    from poltype.config import PoltypeConfig
    from poltype.qm.factory import select_backend

    cfg = PoltypeConfig()
    backend = select_backend(cfg)
    # backend is a Psi4Backend, GaussianBackend, or PySCFBackend instance
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from poltype.config.schema import PoltypeConfig

from poltype.qm.backend import QMBackend


def select_backend(
    config: "PoltypeConfig",
    work_dir: Optional[Path] = None,
) -> QMBackend:
    """Instantiate the correct :class:`QMBackend` based on *config*.

    Selection logic:

    1. If ``config.qm.use_gaus`` or ``config.qm.use_gaus_opt_only`` is
       ``True`` → :class:`~poltype.qm.gaussian_backend.GaussianBackend`
       (legacy flag override).
    2. Otherwise the ``config.qm.backend`` field is used:
       ``"gaussian"`` → :class:`~poltype.qm.gaussian_backend.GaussianBackend`,
       ``"pyscf"`` → :class:`~poltype.qm.pyscf_backend.PySCFBackend`,
       ``"psi4"`` (default) → :class:`~poltype.qm.psi4_backend.Psi4Backend`.

    Parameters
    ----------
    config:
        Fully-populated :class:`~poltype.config.schema.PoltypeConfig`.
    work_dir:
        Optional scratch directory override.  If ``None``, each backend
        will use its own default.

    Returns
    -------
    QMBackend
        An instance of the selected concrete backend.
    """
    from poltype.qm.gaussian_backend import GaussianBackend
    from poltype.qm.pyscf_backend import PySCFBackend
    from poltype.qm.psi4_backend import Psi4Backend

    qm = config.qm
    res = config.resources

    num_proc = res.num_proc or 1
    max_mem_gb = res.max_mem_gb

    # Legacy flags override the ``backend`` field for backward compatibility.
    if qm.use_gaus or qm.use_gaus_opt_only:
        return GaussianBackend(
            work_dir=work_dir,
            num_proc=num_proc,
            max_mem_gb=max_mem_gb,
        )

    # Use the explicit ``backend`` field as the primary selector.
    if qm.backend == "gaussian":
        return GaussianBackend(
            work_dir=work_dir,
            num_proc=num_proc,
            max_mem_gb=max_mem_gb,
        )

    if qm.backend == "pyscf":
        return PySCFBackend(
            work_dir=work_dir,
            num_proc=num_proc,
            max_mem_gb=max_mem_gb,
        )

    # Default: Psi4 (backend == "psi4")
    return Psi4Backend(
        work_dir=work_dir,
        psi4_args=qm.psi4_args,
        num_proc=num_proc,
        max_mem_gb=max_mem_gb,
    )
