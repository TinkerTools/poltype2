"""
tests/unit/test_config.py – unit tests for poltype.config.

These tests verify:
- PoltypeConfig can be instantiated with defaults.
- load_config() correctly parses a poltype.ini string.
- Typed coercions (int, float, bool, list) work correctly.
- Unknown keys are silently ignored.
"""

from __future__ import annotations

import os
import textwrap
from pathlib import Path

import pytest

from poltype.config.schema import PoltypeConfig, QMConfig, ResourceConfig
from poltype.config.loader import load_config, _parse_ini_lines


# ---------------------------------------------------------------------------
# Schema tests
# ---------------------------------------------------------------------------


class TestPoltypeConfigDefaults:
    """PoltypeConfig should have sensible defaults without any arguments."""

    def test_instantiation(self):
        cfg = PoltypeConfig()
        assert cfg is not None

    def test_qm_sub_config_type(self):
        cfg = PoltypeConfig()
        assert isinstance(cfg.qm, QMConfig)

    def test_resource_sub_config_type(self):
        cfg = PoltypeConfig()
        assert isinstance(cfg.resources, ResourceConfig)

    def test_qm_defaults(self):
        qm = QMConfig()
        assert qm.opt_method == "MP2"
        assert qm.opt_basis_set == "6-31G*"
        assert qm.backend == "psi4"
        assert qm.use_gaus is False
        assert qm.dont_use_pyscf is False

    def test_resource_defaults(self):
        res = ResourceConfig()
        assert res.num_proc is None
        assert res.max_mem_gb == 16
        assert res.consumption_ratio == 0.8

    def test_force_field_default(self):
        cfg = PoltypeConfig()
        assert cfg.force_field == "AMOEBA"

    def test_prm_file_path_is_absolute(self):
        cfg = PoltypeConfig()
        assert os.path.isabs(cfg.prm_file_path)

    def test_structure_default_is_none(self):
        cfg = PoltypeConfig()
        assert cfg.structure is None


# ---------------------------------------------------------------------------
# _parse_ini_lines tests
# ---------------------------------------------------------------------------


class TestParseIniLines:
    """Low-level ini line parser."""

    def test_simple_assignment(self):
        lines = ["optmethod = wB97X-D\n"]
        raw = _parse_ini_lines(lines)
        assert raw["optmethod"] == "wB97X-D"

    def test_comment_line_ignored(self):
        lines = ["# this is a comment\n", "optmethod = MP2\n"]
        raw = _parse_ini_lines(lines)
        assert "optmethod" in raw
        assert len(raw) == 1

    def test_inline_comment_stripped(self):
        lines = ["optmethod = MP2 # some comment\n"]
        raw = _parse_ini_lines(lines)
        assert raw["optmethod"] == "MP2"

    def test_bare_keyword_is_none(self):
        lines = ["debugmode\n"]
        raw = _parse_ini_lines(lines)
        assert "debugmode" in raw
        assert raw["debugmode"] is None

    def test_none_value_excluded(self):
        lines = ["optmethod = None\n"]
        raw = _parse_ini_lines(lines)
        assert "optmethod" not in raw

    def test_empty_line_ignored(self):
        lines = ["\n", "  \n", "optmethod = MP2\n"]
        raw = _parse_ini_lines(lines)
        assert len(raw) == 1

    def test_keys_lowercased(self):
        lines = ["OptMethod = MP2\n"]
        raw = _parse_ini_lines(lines)
        assert "optmethod" in raw


# ---------------------------------------------------------------------------
# load_config integration tests
# ---------------------------------------------------------------------------


class TestLoadConfig:
    """load_config() should populate PoltypeConfig from a file."""

    def _write_ini(self, tmp_path: Path, content: str) -> Path:
        ini = tmp_path / "poltype.ini"
        ini.write_text(textwrap.dedent(content))
        return ini

    def test_minimal_config(self, tmp_path):
        ini = self._write_ini(tmp_path, "")
        cfg = load_config(ini)
        assert isinstance(cfg, PoltypeConfig)

    def test_opt_method_parsed(self, tmp_path):
        ini = self._write_ini(tmp_path, "optmethod = wB97X-D\n")
        cfg = load_config(ini)
        assert cfg.qm.opt_method == "wB97X-D"

    def test_numproc_parsed(self, tmp_path):
        ini = self._write_ini(tmp_path, "numproc = 4\n")
        cfg = load_config(ini)
        assert cfg.resources.num_proc == 4

    def test_totalcharge_parsed(self, tmp_path):
        ini = self._write_ini(tmp_path, "totalcharge = -1\n")
        cfg = load_config(ini)
        assert cfg.total_charge == -1

    def test_debugmode_bare_keyword(self, tmp_path):
        ini = self._write_ini(tmp_path, "debugmode\n")
        cfg = load_config(ini)
        assert cfg.debug_mode is True

    def test_use_gaus_flag(self, tmp_path):
        ini = self._write_ini(tmp_path, "use_gaus = True\n")
        cfg = load_config(ini)
        assert cfg.qm.use_gaus is True

    def test_vdwprmtypestofit_list(self, tmp_path):
        ini = self._write_ini(tmp_path, "vdwprmtypestofit = S, T, R\n")
        cfg = load_config(ini)
        assert cfg.vdw_prm_types_to_fit == ["S", "T", "R"]

    def test_onlyvdwatomlist_int_list(self, tmp_path):
        ini = self._write_ini(tmp_path, "onlyvdwatomlist = 1, 3, 5\n")
        cfg = load_config(ini)
        assert cfg.only_vdw_atom_list == [1, 3, 5]

    def test_unknown_key_ignored(self, tmp_path):
        ini = self._write_ini(
            tmp_path, "totally_unknown_key = something_weird\n"
        )
        # Should not raise
        cfg = load_config(ini)
        assert isinstance(cfg, PoltypeConfig)

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_config("/nonexistent/path/poltype.ini")

    def test_comment_and_assignment_mixed(self, tmp_path):
        ini = self._write_ini(
            tmp_path,
            """\
            # QM settings
            optmethod = wB97X-D
            # resource settings
            numproc = 8
            """,
        )
        cfg = load_config(ini)
        assert cfg.qm.opt_method == "wB97X-D"
        assert cfg.resources.num_proc == 8

    def test_optloose_sets_convergence(self, tmp_path):
        ini = self._write_ini(tmp_path, "optloose = True\n")
        cfg = load_config(ini)
        assert cfg.qm.opt_loose is True
        assert cfg.qm.opt_convergence == "LOOSE"
