"""
poltype.config.loader – parse ``poltype.ini`` into :class:`PoltypeConfig`.

The ``load_config`` function reads a ``poltype.ini`` file (the same
format understood by the legacy ``PolarizableTyper.__post_init__``
parser) and returns a fully-typed :class:`~poltype.config.schema.PoltypeConfig`.

This replaces the ~400-line ``elif`` chain in the legacy code.  The
implementation is intentionally simple: it builds a flat key→value dict
from the ini file, then uses a registry of typed coercion functions to
populate the dataclass fields.

All field names follow the ``PoltypeConfig`` naming convention (snake_case)
and are mapped from the original camelCase / mixed-case keyword names
used in ``poltype.ini`` files.

Usage::

    from poltype.config import load_config

    cfg = load_config("poltype.ini")
    print(cfg.qm.opt_method)
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from poltype.config.schema import PoltypeConfig, QMConfig, ResourceConfig


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------


def _parse_bool(value: str) -> bool:
    """Parse a string representation of a boolean."""
    return value.strip().lower() not in ("false", "0", "no", "none", "")


def _parse_int(value: str) -> int:
    return int(value.strip())


def _parse_float(value: str) -> float:
    return float(value.strip())


def _parse_str(value: str) -> str:
    return value.strip()


def _parse_str_list(value: str) -> List[str]:
    return [s.strip() for s in value.split(",") if s.strip()]


def _parse_int_list(value: str) -> List[int]:
    return [int(s.strip()) for s in value.split(",") if s.strip()]


# ---------------------------------------------------------------------------
# Raw ini parsing
# ---------------------------------------------------------------------------


def _parse_ini_lines(lines: List[str]) -> Dict[str, Optional[str]]:
    """Return a mapping of lowercase keyword → raw value string.

    Handles:
    - comment lines starting with ``#``
    - inline comments after a ``#``
    - ``keyword = value`` and bare ``keyword`` (flag-only) lines
    """
    result: Dict[str, Optional[str]] = {}
    for raw in lines:
        parts = raw.split()
        if not parts or parts[0].startswith("#"):
            continue
        # strip inline comment
        for i, tok in enumerate(parts):
            if tok.startswith("#"):
                parts = parts[:i]
                break
        if not parts:
            continue
        line = " ".join(parts)
        if "=" in line:
            key_part, _, val_part = line.partition("=")
            key = key_part.strip().lower()
            val = val_part.strip()
            if val.lower() == "none":
                continue  # explicit None → keep default
            result[key] = val
        else:
            # bare keyword → boolean flag
            result[line.strip().lower()] = None
    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_config(
    ini_path: str | Path = "poltype.ini",
    *,
    defaults: Optional[PoltypeConfig] = None,
) -> PoltypeConfig:
    """Parse a ``poltype.ini`` file and return a :class:`PoltypeConfig`.

    Parameters
    ----------
    ini_path:
        Path to the ``poltype.ini`` file.  Defaults to ``poltype.ini`` in
        the current working directory (legacy behaviour).
    defaults:
        Optional base ``PoltypeConfig`` to use as starting point.  If not
        provided, a fresh default-valued config is used.

    Returns
    -------
    PoltypeConfig
        Fully populated, typed configuration object.
    """
    ini_path = Path(ini_path)
    if not ini_path.exists():
        raise FileNotFoundError(f"Config file not found: {ini_path}")

    with ini_path.open() as fh:
        lines = fh.readlines()

    raw = _parse_ini_lines(lines)
    cfg = defaults if defaults is not None else PoltypeConfig()
    qm = cfg.qm
    res = cfg.resources

    # Build mutable copies of the sub-configs so we can patch them.
    # (They are regular dataclasses, so direct attribute assignment works.)

    def _bool(k: str) -> bool:
        v = raw.get(k)
        if v is None:
            return True  # bare keyword → True
        return _parse_bool(v)

    def _get(k: str, default: Any) -> Any:
        return raw[k] if k in raw else default

    # ------------------------------------------------------------------
    # Structure / input file
    # ------------------------------------------------------------------
    if "structure" in raw:
        cfg.structure = Path(_get("structure", cfg.structure))

    # ------------------------------------------------------------------
    # QM settings
    # ------------------------------------------------------------------
    if "optmethod" in raw:
        qm.opt_method = _parse_str(raw["optmethod"])
    if "optbasisset" in raw:
        qm.opt_basis_set = _parse_str(raw["optbasisset"])
    if "optbasissetfile" in raw:
        qm.opt_basis_set_file = _parse_str(raw["optbasissetfile"])
    if "optconvergence" in raw:
        qm.opt_convergence = _parse_str(raw["optconvergence"]).upper()
    if "optloose" in raw:
        qm.opt_loose = _bool("optloose")
        if qm.opt_loose:
            qm.opt_convergence = "LOOSE"
    if "dmabasissetfile" in raw:
        qm.dma_basis_set_file = _parse_str(raw["dmabasissetfile"])
    if "espmethod" in raw:
        qm.esp_method = _parse_str(raw["espmethod"])
    if "espbasisset" in raw:
        qm.esp_basis_set = _parse_str(raw["espbasisset"])
    if "espbasissetfile" in raw:
        qm.esp_basis_set_file = _parse_str(raw["espbasissetfile"])
    if "toroptmethod" in raw:
        qm.tor_opt_method = _parse_str(raw["toroptmethod"])
    if "toroptbasisset" in raw:
        qm.tor_opt_basis_set = _parse_str(raw["toroptbasisset"])
    if "toroptbasissetfile" in raw:
        qm.tor_opt_basis_set_file = _parse_str(raw["toroptbasissetfile"])
    if "torspbasisset" in raw:
        qm.tor_sp_basis_set = _parse_str(raw["torspbasisset"])
    if "torspbasissetfile" in raw:
        qm.tor_sp_basis_set_file = _parse_str(raw["torspbasissetfile"])
    if "torspbasissethalogen" in raw:
        # legacy name maps to the tor_sp_basis_set for halogens
        qm.iodine_tor_sp_basis_set = _parse_str(raw["torspbasissethalogen"])
    if "iodineoptbasisset" in raw:
        qm.iodine_opt_basis_set = _parse_str(raw["iodineoptbasisset"])
    if "iodineoptbasissetfile" in raw:
        qm.iodine_opt_basis_set_file = _parse_str(raw["iodineoptbasissetfile"])
    if "iodinedmabasissetfile" in raw:
        qm.iodine_dma_basis_set_file = _parse_str(raw["iodinedmabasissetfile"])
    if "iodineespbasissetfile" in raw:
        qm.iodine_esp_basis_set_file = _parse_str(raw["iodineespbasissetfile"])
    if "iodinetoroptbasisset" in raw:
        qm.iodine_tor_opt_basis_set = _parse_str(raw["iodinetoroptbasisset"])
    if "iodinetoroptbasissetfile" in raw:
        qm.iodine_tor_opt_basis_set_file = _parse_str(
            raw["iodinetoroptbasissetfile"]
        )
    if "iodinetorspbasisset" in raw:
        qm.iodine_tor_sp_basis_set = _parse_str(raw["iodinetorspbasisset"])
    if "iodinetorspbasissetfile" in raw:
        qm.iodine_tor_sp_basis_set_file = _parse_str(
            raw["iodinetorspbasissetfile"]
        )
    if "use_gaus" in raw and "opt" not in raw.get("use_gaus", ""):
        qm.use_gaus = _bool("use_gaus")
    if "use_gausoptonly" in raw:
        qm.use_gaus_opt_only = _bool("use_gausoptonly")
    if "use_psi4_geometric_opt" in raw:
        qm.use_psi4_geometric_opt = _bool("use_psi4_geometric_opt")
    if "dont_use_pyscf" in raw:
        qm.dont_use_pyscf = _bool("dont_use_pyscf")
    if "psi4_args" in raw:
        qm.psi4_args = _parse_str(raw["psi4_args"])
    if "toroptmethodlist" in raw:
        qm.tor_opt_method_list = _parse_str_list(raw["toroptmethodlist"])
    if "torspmethodlist" in raw:
        qm.tor_sp_method_list = _parse_str_list(raw["torspmethodlist"])
    if "sameleveldmaesp" in raw:
        qm.same_level_dma_esp = _bool("sameleveldmaesp")

    # ------------------------------------------------------------------
    # Resource settings
    # ------------------------------------------------------------------
    if "numproc" in raw:
        res.num_proc = _parse_int(raw["numproc"])
    if "maxmem" in raw:
        # accept e.g. "16GB" or "16"
        val = raw["maxmem"].replace("GB", "").replace("gb", "").strip()
        res.max_mem_gb = int(val)
    if "maxdisk" in raw:
        val = raw["maxdisk"].replace("GB", "").replace("gb", "").strip()
        res.max_disk_gb = int(val)
    if "jobsatsametime" in raw and "maxjobsatsametime" not in raw and "parentjobsatsametime" not in raw:
        res.jobs_at_same_time = _parse_int(raw["jobsatsametime"])
    if "maxjobsatsametime" in raw:
        res.max_jobs_at_same_time = _parse_float(raw["maxjobsatsametime"])
    if "parentjobsatsametime" in raw:
        res.parent_jobs_at_same_time = _parse_int(raw["parentjobsatsametime"])
    if "coresperjob" in raw:
        res.cores_per_job = _parse_int(raw["coresperjob"])
    if "maximizejobsatsametime" in raw:
        res.maximize_jobs_at_same_time = _bool("maximizejobsatsametime")
    if "consumptionratio" in raw:
        res.consumption_ratio = _parse_float(raw["consumptionratio"])

    # ------------------------------------------------------------------
    # Force field
    # ------------------------------------------------------------------
    if "prmfilepath" in raw:
        cfg.prm_file_path = _parse_str(raw["prmfilepath"])

    # ------------------------------------------------------------------
    # Pipeline control
    # ------------------------------------------------------------------
    if "databasematchonly" in raw:
        cfg.database_match_only = _bool("databasematchonly")
    if "setupfragjobsonly" in raw:
        cfg.setup_frag_jobs_only = _bool("setupfragjobsonly")
    if "optonly" in raw:
        cfg.opt_only = _bool("optonly")
    if "checkinputonly" in raw:
        cfg.check_input_only = _bool("checkinputonly")
    if "generateinputfilesonly" in raw:
        cfg.generate_input_files_only = _bool("generateinputfilesonly")
    if "genprotstatesonly" in raw:
        cfg.gen_prot_states_only = _bool("genprotstatesonly")

    # ------------------------------------------------------------------
    # Molecule / protonation
    # ------------------------------------------------------------------
    if "totalcharge" in raw:
        cfg.total_charge = _parse_int(raw["totalcharge"])
    if "addhydrogens" in raw:
        cfg.add_hydrogens = _bool("addhydrogens")
    if "addhydrogentocharged" in raw:
        cfg.add_hydrogen_to_charged = _bool("addhydrogentocharged")
    if "allowradicals" in raw:
        cfg.allow_radicals = _bool("allowradicals")

    # ------------------------------------------------------------------
    # Conformer generation
    # ------------------------------------------------------------------
    if "generateextendedconf" in raw:
        cfg.generate_extended_conf = _bool("generateextendedconf")
    if "generate_symm_frag_conf" in raw:
        cfg.generate_symm_frag_conf = _bool("generate_symm_frag_conf")
    if "userconformation" in raw:
        cfg.user_conformation = _bool("userconformation")
    if "numespconfs" in raw:
        cfg.num_esp_confs = _parse_int(raw["numespconfs"])
    if "espextraconflist" in raw:
        cfg.esp_extra_conf_list = _parse_str_list(raw["espextraconflist"])

    # ------------------------------------------------------------------
    # Symmetry
    # ------------------------------------------------------------------
    if "usesymtypes" in raw:
        cfg.use_sym_types = _bool("usesymtypes")

    # ------------------------------------------------------------------
    # Torsion settings
    # ------------------------------------------------------------------
    if "tortor" in raw:
        cfg.tor_tor = _bool("tortor")
    if "torfit2drotonly" in raw:
        cfg.tor_fit_2d_rot_only = _bool("torfit2drotonly")
    if "torfit1drotonly" in raw:
        cfg.tor_fit_1d_rot_only = _bool("torfit1drotonly")
    if "defaultmaxtorsiongridpoints" in raw:
        cfg.default_max_torsion_grid_points = _parse_int(
            raw["defaultmaxtorsiongridpoints"]
        )
    if "maxtorrmspdrrel" in raw:
        cfg.max_tor_rms_pdr_rel = _parse_float(raw["maxtorrmspdrrel"])
    if "rotalltors" in raw:
        cfg.rot_all_tors = _bool("rotalltors")
    if "onlyrottortorlist" in raw:
        cfg.only_rot_tor_tor_list = _parse_str_list(raw["onlyrottortorlist"])
    if "dontrotbndslist" in raw:
        cfg.dont_rot_bnds_list = _parse_str_list(raw["dontrotbndslist"])
    if "nonaroringtor1dscan" in raw:
        cfg.non_aro_ring_tor_1d_scan = _bool("nonaroringtor1dscan")
    if "maxtorresnitrogen" in raw:
        cfg.max_tor_res_nitrogen = _parse_int(raw["maxtorresnitrogen"])
    if "xtbtorresconstant" in raw:
        cfg.xtb_tor_res_constant = _parse_float(raw["xtbtorresconstant"])

    # ------------------------------------------------------------------
    # ESP / multipole
    # ------------------------------------------------------------------
    if "esprestweight" in raw:
        cfg.esp_rest_weight = _parse_float(raw["esprestweight"])
    if "espgrad" in raw:
        cfg.esp_grad = _parse_float(raw["espgrad"])
    if "potentialoffset" in raw:
        cfg.potential_offset = _parse_float(raw["potentialoffset"])
    if "skipespfiterror" in raw:
        cfg.skip_esp_fit_error = _bool("skipespfiterror")
    if "fitfirsttorsionfoldphase" in raw:
        cfg.fit_first_torsion_fold_phase = _bool("fitfirsttorsionfoldphase")
    if "skipgridsearch" in raw:
        cfg.skip_grid_search = _bool("skipgridsearch")

    # ------------------------------------------------------------------
    # vdW
    # ------------------------------------------------------------------
    if "dovdwscan" in raw:
        cfg.do_vdw_scan = _bool("dovdwscan")
    if "vdwprmtypestofit" in raw:
        cfg.vdw_prm_types_to_fit = _parse_str_list(raw["vdwprmtypestofit"])
    if "fixvdwtyperadii" in raw:
        cfg.fix_vdw_type_radii = _parse_str_list(raw["fixvdwtyperadii"])
    if "addlonepairvdwsites" in raw:
        cfg.add_lone_pair_vdw_sites = _bool("addlonepairvdwsites")
    if "accuratevdwsp" in raw:
        cfg.accurate_vdw_sp = _bool("accuratevdwsp")
    if "use_qmopt_vdw" in raw:
        cfg.use_qmopt_vdw = _bool("use_qmopt_vdw")
    if "use_gau_vdw" in raw:
        cfg.use_gau_vdw = _bool("use_gau_vdw")
    if "vdwmaxqmstartingpointspertype" in raw:
        cfg.vdw_max_qm_starting_points_per_type = _parse_int(
            raw["vdwmaxqmstartingpointspertype"]
        )
    if "vdwmaxtinkergridpoints" in raw:
        cfg.vdw_max_tinker_grid_points = _parse_int(
            raw["vdwmaxtinkergridpoints"]
        )
    if "homodimers" in raw:
        cfg.homo_dimers = _bool("homodimers")
    if "onlyvdwatomlist" in raw:
        cfg.only_vdw_atom_list = _parse_int_list(raw["onlyvdwatomlist"])

    # ------------------------------------------------------------------
    # Fragmentation
    # ------------------------------------------------------------------
    if "smallmoleculefragmenter" in raw:
        cfg.small_molecule_fragmenter = _bool("smallmoleculefragmenter")
    if "quickdatabasesearch" in raw:
        cfg.quick_database_search = _bool("quickdatabasesearch")
    if "fitred" in raw:
        cfg.fit_red = _bool("fitred")

    # ------------------------------------------------------------------
    # Output / debug
    # ------------------------------------------------------------------
    if "writeoutmultipole" in raw:
        cfg.write_out_multipole = _bool("writeoutmultipole")
    if "writeoutbond" in raw:
        cfg.write_out_bond = _bool("writeoutbond")
    if "writeoutangle" in raw:
        cfg.write_out_angle = _bool("writeoutangle")
    if "writeoutpolarize" in raw:
        cfg.write_out_polarize = _bool("writeoutpolarize")
    if "writeouttorsion" in raw:
        cfg.write_out_torsion = _bool("writeouttorsion")
    if "debugmode" in raw:
        cfg.debug_mode = _bool("debugmode")
    if "tordebugmode" in raw:
        cfg.tor_debug_mode = _bool("tordebugmode")
    if "fragmenterdebugmode" in raw:
        cfg.fragmenter_debug_mode = _bool("fragmenterdebugmode")
    if "toroptdebugmode" in raw:
        cfg.tor_opt_debug_mode = _bool("toroptdebugmode")
    if "deleteallnonqmfiles" in raw:
        cfg.delete_all_non_qm_files = _bool("deleteallnonqmfiles")
    if "lastlogfileupdatetime" in raw:
        cfg.last_log_file_update_time = _parse_float(
            raw["lastlogfileupdatetime"]
        )
    if "useuniquefilenames" in raw:
        cfg.use_unique_filenames = _bool("useuniquefilenames")

    # ------------------------------------------------------------------
    # ML potentials
    # ------------------------------------------------------------------
    if "anienvname" in raw:
        cfg.ani_env_name = _parse_str(raw["anienvname"])
    if "xtbenvname" in raw:
        cfg.xtb_env_name = _parse_str(raw["xtbenvname"])
    if "anifmax" in raw:
        cfg.ani_fmax = _parse_float(raw["anifmax"])
    if "xtbmethod" in raw:
        cfg.xtb_method = _parse_int(raw["xtbmethod"])
    if "fennixmodelname" in raw:
        cfg.fennix_model_name = _parse_str(raw["fennixmodelname"])

    # ------------------------------------------------------------------
    # Misc
    # ------------------------------------------------------------------
    if "polareps" in raw:
        cfg.polar_eps = _parse_str(raw["polareps"])
    if "boltzmantemp" in raw:
        cfg.boltzman_temp = _parse_float(raw["boltzmantemp"])
    if "absdipoletol" in raw:
        cfg.abs_dipole_tol = _parse_float(raw["absdipoletol"])
    if "mmbondtol" in raw:
        cfg.mm_bond_tol = _parse_float(raw["mmbondtol"])
    if "mmangletol" in raw:
        cfg.mm_angle_tol = _parse_float(raw["mmangletol"])
    if "transferanyhydrogentor" in raw:
        cfg.transfer_any_hydrogen_tor = _bool("transferanyhydrogentor")
    if "pcm_auto" in raw:
        cfg.pcm_auto = _bool("pcm_auto")
    if "usepoleditframes" in raw:
        cfg.use_poledit_frames = _bool("usepoleditframes")
    if "indextotypefile" in raw:
        cfg.index_to_type_file = _parse_str(raw["indextotypefile"])
    if "indextompoleframefile" in raw:
        cfg.index_to_mpole_frame_file = _parse_str(
            raw["indextompoleframefile"]
        )
    if "bashrcpath" in raw:
        cfg.bashrc_path = _parse_str(raw["bashrcpath"])
    if "parentname" in raw:
        cfg.parent_name = _parse_str(raw["parentname"])
    if "checktraj" in raw:
        cfg.check_traj = _bool("checktraj")
    if "inputkeyfile" in raw:
        src = _parse_str(raw["inputkeyfile"])
        safe = "inputkey.key"
        if os.path.exists(src):
            shutil.copy(src, safe)
        cfg.input_key_file = safe

    # Write back mutated sub-configs (they are mutable dataclass instances
    # but reassignment makes intent clearer when using frozen=False).
    cfg.qm = qm
    cfg.resources = res

    return cfg
