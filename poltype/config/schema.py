"""
poltype.config.schema – typed configuration dataclasses.

All user-facing Poltype settings live here, grouped into logical
sub-configs.  Each field has an explicit Python type and a documented
default so users and tooling can discover available knobs without
reading the source.

These dataclasses are *immutable by convention*: the pipeline creates
one ``PoltypeConfig`` at startup and never mutates it during a run.
Mutable per-run state belongs in ``poltype.pipeline.context``.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Literal, Optional

# ---------------------------------------------------------------------------
# Helper: resolve a path relative to the ParameterFiles data directory.
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parents[2]
_PARAM_DIR = _REPO_ROOT / "ParameterFiles"
_MODULE_DIR = _REPO_ROOT / "PoltypeModules"


def _param(filename: str) -> str:
    """Return absolute path string for a file under ParameterFiles/."""
    return str(_PARAM_DIR / filename)


def _module(subpath: str) -> str:
    """Return absolute path string for a file under PoltypeModules/."""
    return str(_MODULE_DIR / subpath)


# ---------------------------------------------------------------------------
# QM method / basis-set configuration
# ---------------------------------------------------------------------------


@dataclass
class QMConfig:
    """All quantum-mechanics method and basis-set choices.

    These settings map directly to the keyword options previously parsed
    by the 400-line ``elif`` chain in ``PolarizableTyper.__post_init__``.
    """

    # --- geometry optimisation ---
    opt_method: str = "MP2"
    opt_basis_set: str = "6-31G*"
    opt_basis_set_file: str = "6-31g_st_.0.gbs"
    opt_max_cycles: int = 400
    opt_convergence: str = "LOOSE"
    opt_loose: bool = True

    # --- DMA multipole ---
    dma_basis_set: str = "6-311G**"
    dma_basis_set_file: str = "6-311g_st__st_.0.gbs"

    # --- ESP fitting ---
    esp_method: str = "MP2"
    esp_basis_set: str = "aug-cc-pVTZ"
    esp_basis_set_file: str = "aug-cc-pvtz.1.gbs"

    # --- torsion optimisation ---
    tor_opt_method: str = "wB97X-D"
    tor_opt_basis_set: str = "6-31G*"
    tor_opt_basis_set_file: str = "6-31g_st_.0.gbs"

    # --- torsion single-point ---
    tor_sp_method: str = "wB97X-D"
    tor_sp_basis_set: str = "6-311+G*"
    tor_sp_basis_set_file: str = "6-311+g_st_.0.gbs"

    # --- halogen / iodine overrides ---
    iodine_opt_basis_set: str = "def2-svp"
    iodine_opt_basis_set_file: str = "def2-svp.1.gbs"
    iodine_dma_basis_set_file: str = "def2-svp.1.gbs"
    iodine_esp_basis_set_file: str = "def2-tzvpp.1.gbs"
    iodine_tor_opt_basis_set: str = "def2-svp"
    iodine_tor_opt_basis_set_file: str = "def2-svp.1.gbs"
    iodine_tor_sp_basis_set: str = "def2-svp"
    iodine_tor_sp_basis_set_file: str = "def2-svp.1.gbs"

    # --- backend selection ---
    backend: Literal["psi4", "gaussian", "pyscf"] = "psi4"
    use_gaus: bool = False
    use_gaus_opt_only: bool = False
    use_psi4_geometric_opt: bool = True
    dont_use_pyscf: bool = False

    # --- Psi4-specific ---
    psi4_args: str = " --loglevel 30 "
    same_level_dma_esp: bool = False

    # --- method-list overrides (per-torsion-job) ---
    tor_opt_method_list: List[str] = field(default_factory=list)
    tor_sp_method_list: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Compute resource configuration
# ---------------------------------------------------------------------------


@dataclass
class ResourceConfig:
    """CPU / memory / scratch limits for the run."""

    num_proc: Optional[int] = None  # ``None`` → auto-detect
    max_mem_gb: int = 16
    max_disk_gb: int = 100
    scratch_dir: Path = Path("/tmp")

    # parallelism knobs
    jobs_at_same_time: int = 0  # 0 → auto
    max_jobs_at_same_time: float = 10
    parent_jobs_at_same_time: int = 1
    cores_per_job: int = 2
    maximize_jobs_at_same_time: bool = True
    consumption_ratio: float = 0.8


# ---------------------------------------------------------------------------
# Top-level user configuration
# ---------------------------------------------------------------------------


@dataclass
class PoltypeConfig:
    """All user-facing Poltype settings in one typed, validated object.

    This replaces the ~400-line ``elif`` keyword parser in
    ``PolarizableTyper.__post_init__``.  The ``load_config`` function in
    ``poltype.config.loader`` populates this from a ``poltype.ini`` file.

    Only a ``structure`` path is required; everything else has a sensible
    default.
    """

    # --- required ---
    structure: Optional[Path] = None  # input molecule file (.sdf/.mol/.mol2)

    # --- sub-configs ---
    qm: QMConfig = field(default_factory=QMConfig)
    resources: ResourceConfig = field(default_factory=ResourceConfig)

    # --- force-field ---
    force_field: Literal["AMOEBA", "AMOEBA+"] = "AMOEBA"
    prm_file_path: str = _param("amoebabio18.prm")
    parameter_files_path: str = str(_PARAM_DIR)

    # AMOEBA+ paths
    amoeba_plus_nonbonded_prm: str = _module(
        "lDatabaseParser/prm/amoebaplusNonbonded.prm"
    )
    amoeba_plus_nonbonded_dat: str = _module(
        "lDatabaseParser/dat/amoebaplusNonbondedType.dat"
    )
    ldatabaseparser_path: str = _module(
        "lDatabaseParser/lAssignAMOEBAplusPRM.py"
    )
    ldatabaseparser_prm_dir: str = _module("lDatabaseParser/prm")

    # AMOEBA small-molecule libraries
    small_molecule_prm_lib: str = _param("amoeba09.prm")
    latest_small_molecule_prm_lib: str = _param("amoeba21.prm")
    latest_small_molecule_polarize_prm_lib: str = _param("amoeba21polarize.prm")
    updated_small_molecule_polarize_prm_lib: str = _module(
        "lDatabaseParser/prm/polarize.prm"
    )
    small_molecule_mm3_prm_lib: str = _param("mm3.prm")

    # SMARTS mapping files
    small_molecule_smarts_to_tinker_class: str = _param(
        "amoeba21smartstoclass.txt"
    )
    small_molecule_smarts_to_tinker_descrip: str = _param(
        "smartstoamoebatypedescrip.txt"
    )
    small_molecule_smarts_to_mm3_descrip: str = _param(
        "smartstomm3typedescrip.txt"
    )
    latest_small_molecule_smarts_to_types_polarize: str = _param(
        "amoeba21polarcommenttoparameters.txt"
    )
    smarts_to_solute_radii_map: str = _param("SMARTsToSoluteRadiiMap.txt")
    external_parameter_database: str = _param("externalparameterdatabase.txt")

    # Tinker / GDMA executables
    analyze_path: str = "analyze"
    poledit_path: str = "poledit"
    potential_path: str = "potential"
    minimize_path: str = "minimize"
    dynamic_path: str = "dynamic"
    bar_path: str = "bar"
    gdma_path: str = "gdma"

    # --- pipeline control flags ---
    do_torsion_fit: bool = True
    tor_fit: bool = True
    database_match_only: bool = False
    setup_frag_jobs_only: bool = False
    opt_only: bool = False
    check_input_only: bool = False
    generate_input_files_only: bool = False
    gen_prot_states_only: bool = False
    dry_run: bool = False
    skip_validation: bool = False

    # --- molecule / protonation ---
    total_charge: Optional[int] = None
    add_hydrogens: bool = False
    add_hydrogen_to_charged: bool = True
    allow_radicals: bool = False

    # --- conformer generation ---
    generate_extended_conf: bool = True
    generate_symm_frag_conf: bool = False
    user_conformation: bool = False
    num_esp_confs: int = 1
    esp_extra_conf_list: List[str] = field(default_factory=list)

    # --- symmetry ---
    use_sym_types: bool = True

    # --- torsion settings ---
    tor_tor: bool = False
    tor_fit_2d_rot_only: bool = False
    tor_fit_1d_rot_only: bool = False
    fold_num: int = 6
    default_max_torsion_grid_points: int = 40
    max_tor_rms_pdr_rel: float = 0.1
    rot_all_tors: bool = False
    only_rot_tor_tor_list: List[str] = field(default_factory=list)
    dont_rot_bnds_list: List[str] = field(default_factory=list)
    non_aro_ring_tor_1d_scan: bool = False
    max_tor_res_nitrogen: int = 2
    xtb_tor_res_constant: float = 5.0

    # --- ESP / multipole settings ---
    esp_rest_weight: float = 1.0
    esp_grad: float = 0.1
    potential_offset: float = 1.0
    skip_esp_fit_error: bool = False
    fit_first_torsion_fold_phase: bool = False
    skip_grid_search: bool = True

    # --- vdW settings ---
    do_vdw_scan: bool = False
    vdw_prm_types_to_fit: List[str] = field(default_factory=lambda: ["S", "T"])
    fix_vdw_type_radii: List[str] = field(default_factory=list)
    add_lone_pair_vdw_sites: bool = False
    accurate_vdw_sp: bool = False
    use_qmopt_vdw: bool = False
    use_gau_vdw: bool = False
    vdw_max_qm_starting_points_per_type: int = 1
    vdw_max_tinker_grid_points: int = 50
    homo_dimers: bool = False
    only_vdw_atom_list: Optional[List[int]] = None

    # --- fragmentation ---
    use_fragmentation: bool = True
    small_molecule_fragmenter: bool = False
    wbo_tolerance: float = 0.05
    max_growth_cycles: int = 5
    quick_database_search: bool = False
    fit_red: bool = False

    # --- output / debug ---
    write_out_multipole: bool = True
    write_out_bond: bool = True
    write_out_angle: bool = True
    write_out_polarize: bool = True
    write_out_torsion: bool = True
    debug_mode: bool = False
    tor_debug_mode: bool = False
    fragmenter_debug_mode: bool = False
    tor_opt_debug_mode: bool = False
    delete_all_non_qm_files: bool = False
    deleted_files: bool = False
    last_log_file_update_time: float = 1.0
    use_unique_filenames: bool = False

    # --- ML potentials ---
    ani_env_name: str = "ani"
    xtb_env_name: str = "xtbenv"
    ani_fmax: float = 0.05
    xtb_method: int = 2
    fennix_model_dir: str = str(_REPO_ROOT / "FennixModels")
    fennix_model_name: str = "ani2x_model0"

    # --- misc ---
    polar_eps: Optional[str] = None
    boltzman_temp: float = 8.0
    abs_dipole_tol: float = 0.5
    mm_bond_tol: float = 0.5
    mm_angle_tol: float = 0.5
    transfer_any_hydrogen_tor: bool = True
    pcm_auto: bool = True
    use_poledit_frames: bool = True
    input_key_file: Optional[str] = None
    index_to_type_file: Optional[str] = None
    index_to_mpole_frame_file: Optional[str] = None
    polt_type_ini: bool = True  # whether to read poltype.ini from CWD
    bashrc_path: Optional[str] = None
    parent_name: Optional[str] = None
    prm_mod_list: List[str] = field(default_factory=list)
    only_fit_tors_together: List[str] = field(default_factory=list)
    check_traj: bool = False
    opt_method_list: List[str] = field(default_factory=list)
    same_level_dma_esp: bool = False
    num_esp_confs: int = 1
