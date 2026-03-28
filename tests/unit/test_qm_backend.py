"""
tests/unit/test_qm_backend.py – unit tests for poltype.qm abstractions.

Tests cover:
- QMBackend is an ABC that cannot be instantiated directly.
- Concrete stub backends (Psi4, Gaussian, PySCF) satisfy the interface.
- Stub methods raise NotImplementedError.
- Result dataclasses can be instantiated.
- select_backend() returns the correct backend type based on config.
"""

from __future__ import annotations

import numpy as np
import pytest

from poltype.config.schema import PoltypeConfig, QMConfig, ResourceConfig
from poltype.molecule.molecule import Molecule
from poltype.qm.backend import (
    ESPGridResult,
    OptimizationResult,
    QMBackend,
    TorsionScanResult,
    WBOResult,
)
from poltype.qm.factory import select_backend
from poltype.qm.gaussian_backend import GaussianBackend
from poltype.qm.psi4_backend import Psi4Backend
from poltype.qm.pyscf_backend import PySCFBackend


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def simple_molecule():
    return Molecule.from_smiles("C", name="methane", embed_3d=True)


@pytest.fixture
def default_config():
    return PoltypeConfig()


# ---------------------------------------------------------------------------
# ABC enforcement
# ---------------------------------------------------------------------------


class TestQMBackendABC:
    def test_cannot_instantiate_abc(self):
        with pytest.raises(TypeError):
            QMBackend()  # type: ignore[abstract]

    def test_concrete_must_implement_all_methods(self):
        """A class that only partially implements QMBackend should fail."""

        class PartialBackend(QMBackend):
            def optimize_geometry(self, *a, **kw):
                pass
            # Missing: compute_esp_grid, torsion_scan, compute_wbo_matrix

        with pytest.raises(TypeError):
            PartialBackend()


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------


class TestResultDataclasses:
    def test_optimization_result(self):
        coords = np.zeros((5, 3))
        r = OptimizationResult(coordinates=coords, energy=-100.0)
        assert r.converged is True
        assert r.log_path is None

    def test_esp_grid_result(self):
        grid = np.zeros((100, 3))
        vals = np.zeros(100)
        r = ESPGridResult(grid_coords=grid, esp_values=vals)
        assert r.log_path is None

    def test_torsion_scan_result(self):
        angles = np.arange(-180, 180, 15, dtype=float)
        energies = np.zeros_like(angles)
        r = TorsionScanResult(angles=angles, energies=energies)
        assert len(r.log_paths) == 0

    def test_wbo_result(self):
        mat = np.eye(5)
        r = WBOResult(wbo_matrix=mat)
        assert r.log_path is None


# ---------------------------------------------------------------------------
# Concrete backend interface conformance
# ---------------------------------------------------------------------------


class TestPsi4BackendInterface:
    def test_instantiation(self):
        b = Psi4Backend()
        assert isinstance(b, QMBackend)

    def test_name(self):
        assert Psi4Backend().name == "psi4"

    def test_optimize_geometry_raises(self, simple_molecule):
        b = Psi4Backend()
        with pytest.raises(NotImplementedError):
            b.optimize_geometry(simple_molecule, "MP2", "6-31G*")

    def test_compute_esp_grid_raises(self, simple_molecule):
        b = Psi4Backend()
        with pytest.raises(NotImplementedError):
            b.compute_esp_grid(simple_molecule, "MP2", "aug-cc-pVTZ")

    def test_torsion_scan_raises(self, simple_molecule):
        b = Psi4Backend()
        with pytest.raises(NotImplementedError):
            b.torsion_scan(simple_molecule, (0, 1), "wB97X-D", "6-31G*")

    def test_compute_wbo_matrix_raises(self, simple_molecule):
        b = Psi4Backend()
        with pytest.raises(NotImplementedError):
            b.compute_wbo_matrix(simple_molecule, "HF", "6-31G*")


class TestGaussianBackendInterface:
    def test_instantiation(self):
        b = GaussianBackend()
        assert isinstance(b, QMBackend)

    def test_name(self):
        assert GaussianBackend().name == "gaussian"

    def test_optimize_geometry_raises(self, simple_molecule):
        with pytest.raises(NotImplementedError):
            GaussianBackend().optimize_geometry(simple_molecule, "MP2", "6-31G*")


class TestPySCFBackendInterface:
    def test_instantiation(self):
        b = PySCFBackend()
        assert isinstance(b, QMBackend)

    def test_name(self):
        assert PySCFBackend().name == "pyscf"

    def test_optimize_geometry_raises(self, simple_molecule):
        with pytest.raises(NotImplementedError):
            PySCFBackend().optimize_geometry(simple_molecule, "MP2", "6-31G*")


# ---------------------------------------------------------------------------
# Factory / select_backend tests
# ---------------------------------------------------------------------------


class TestSelectBackend:
    def test_default_is_pyscf(self, default_config):
        # Default: use_gaus=False, dont_use_pyscf=False → PySCF
        backend = select_backend(default_config)
        assert isinstance(backend, PySCFBackend)

    def test_dont_use_pyscf_selects_psi4(self):
        cfg = PoltypeConfig(qm=QMConfig(dont_use_pyscf=True))
        backend = select_backend(cfg)
        assert isinstance(backend, Psi4Backend)

    def test_use_gaus_selects_gaussian(self):
        cfg = PoltypeConfig(qm=QMConfig(use_gaus=True))
        backend = select_backend(cfg)
        assert isinstance(backend, GaussianBackend)

    def test_use_gaus_opt_only_selects_gaussian(self):
        cfg = PoltypeConfig(qm=QMConfig(use_gaus_opt_only=True))
        backend = select_backend(cfg)
        assert isinstance(backend, GaussianBackend)

    def test_psi4_args_propagated(self):
        cfg = PoltypeConfig(
            qm=QMConfig(dont_use_pyscf=True, psi4_args="--test")
        )
        backend = select_backend(cfg)
        assert isinstance(backend, Psi4Backend)
        assert backend.psi4_args == "--test"

    def test_num_proc_propagated(self):
        cfg = PoltypeConfig(
            qm=QMConfig(dont_use_pyscf=True),
            resources=ResourceConfig(num_proc=8),
        )
        backend = select_backend(cfg)
        assert backend.num_proc == 8

    def test_repr_contains_backend_name(self, default_config):
        backend = select_backend(default_config)
        assert backend.name in repr(backend)
