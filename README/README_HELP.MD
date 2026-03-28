## Advanced Usage

All settings are specified in a `poltype.ini` file or via CLI flags.  The
pipeline reads configuration with `poltype.config.loader.load_config()` and
populates a `PoltypeConfig` dataclass (see `poltype/config/schema.py` for the
full list of typed fields and their defaults).

### Table of Contents

- [CLI Options](#cli-options)
- [Required Inputs](#required-inputs)
- [Resource Control](#resource-control)
- [Protonation State](#protonation-state)
- [Quantum Chemistry Packages](#quantum-chemistry-packages)
- [Atom Typing](#atom-typing)
- [Geometry Optimisation](#geometry-optimisation)
- [Electrostatics (DMA and ESP)](#electrostatics-dma-and-esp)
- [Database Matching](#database-matching)
- [Torsion Parameterisation](#torsion-parameterisation)
- [Fragmentation](#fragmentation)
- [Van der Waals Parameterisation](#van-der-waals-parameterisation)
- [Pipeline Control](#pipeline-control)

#### Notes

* All booleans accept `True` or `False` values (bare keywords without `=` are treated as `True`).
* Keywords are case-insensitive.
* Defaults are chosen to maximise accuracy and efficiency.
* Large molecules may require extra memory for QM jobs.

---

### CLI Options

```
python -m poltype [OPTIONS]

  -c, --config FILE        poltype.ini path (default: poltype.ini)
  -w, --work-dir DIR       Working directory (default: current directory)
  -v, --verbose            Enable debug logging
  -q, --quiet              Suppress output except warnings and errors
  --resume                 Resume from the last checkpoint
  --no-checkpoint          Disable checkpoint saving
  --dry-run                Validate inputs without executing QM jobs
  --version                Show version and exit
```

---

### Required Inputs

`structure` — Input structure file (SDF recommended). Must include bond
orders consistent with the desired total charge. **Required.**

---

### Resource Control

| Keyword | Description | Default |
|---------|-------------|---------|
| `numproc` | Number of processors | auto-detect |
| `maxmem` | Maximum RAM in GB (e.g. `16` or `16GB`) | 16 |
| `maxdisk` | Maximum scratch disk in GB | 100 |
| `consumptionratio` | Fraction of node resources to use | 0.8 |
| `coresperjob` | Cores per parallel QM job | 2 |
| `jobsatsametime` | Maximum concurrent fragment jobs | auto |
| `maxjobsatsametime` | Hard cap on concurrent jobs | 10 |
| `parentjobsatsametime` | Number of Poltype jobs on same host | 1 |
| `maximizejobsatsametime` | Maximise number of concurrent jobs | `True` |

---

### Protonation State

| Keyword | Description | Default |
|---------|-------------|---------|
| `addhydrogens` | Add hydrogens to unprotonated atoms | `False` |
| `addhydrogentocharged` | Add hydrogens to under-valent C/N | `True` |
| `allowradicals` | Allow radical (open-shell) species | `False` |
| `genprotstatesonly` | Output protonation states and exit | `False` |

---

### Quantum Chemistry Packages

| Keyword | Description | Default |
|---------|-------------|---------|
| `use_gaus` | Use Gaussian for all QM calculations | `False` |
| `use_gausoptonly` | Use Gaussian only for geometry optimisation | `False` |
| `dont_use_pyscf` | Disable PySCF fallback | `False` |
| `use_psi4_geometric_opt` | Use Psi4 `geomeTRIC` optimiser | `True` |
| `psi4_args` | Extra command-line arguments passed to Psi4 | `--loglevel 30` |

By default, Psi4 is the primary QM engine.  Set `use_gaus=True` to use
Gaussian for all jobs, or `use_gausoptonly=True` to use Gaussian only for
geometry optimisation.

---

### Atom Typing

| Keyword | Description | Default |
|---------|-------------|---------|
| `usesymtypes` | Use symmetry to reduce and average atom types | `True` |
| `indextotypefile` | Custom atom-index-to-type mapping file | — |
| `indextompoleframefile` | Custom atom-index-to-multipole-frame file | — |

---

### Geometry Optimisation

| Keyword | Description | Default |
|---------|-------------|---------|
| `optmethod` | QM method for geometry optimisation | `MP2` |
| `optbasisset` | Basis set for geometry optimisation | `6-31G*` |
| `optbasissetfile` | Basis set file for geometry optimisation | `6-31g_st_.0.gbs` |
| `optconvergence` | Convergence level (`LOOSE`/`VERYLOOSE`) | `LOOSE` |
| `optloose` | Shorthand for `LOOSE` convergence | `True` |
| `generateextendedconf` | Use extended conformation with restrained torsions | `True` |
| `userconformation` | Use user-supplied conformation (skip 3D generation) | `False` |
| `pcm_auto` | Auto-activate PCM for zwitterions | `True` |

---

### Electrostatics (DMA and ESP)

| Keyword | Description | Default |
|---------|-------------|---------|
| `dmabasissetfile` | Basis set file for DMA single-point | `6-311g_st__st_.0.gbs` |
| `espmethod` | QM method for ESP single-point | `MP2` |
| `espbasisset` | Basis set for ESP single-point | `aug-cc-pVTZ` |
| `espbasissetfile` | Basis set file for ESP single-point | `aug-cc-pvtz.1.gbs` |
| `sameleveldmaesp` | Share a single QM job for DMA and ESP | `False` |
| `numespconfs` | Number of conformations for ESP fitting | 1 |
| `espextraconflist` | Comma-separated extra conformer files for ESP | — |
| `esprestweight` | ESP restraint weight | 1.0 |
| `espgrad` | ESP gradient convergence criterion (kcal/mol/e) | 0.1 |
| `potentialoffset` | Grid offset from vdW surface (Å) | 1.0 |
| `usepoleditframes` | Use POLEDIT for multipole frame detection | `True` |

---

### Database Matching

| Keyword | Description | Default |
|---------|-------------|---------|
| `writeoutmultipole` | Write multipole terms from database | `True` |
| `writeoutbond` | Write bond terms from database | `True` |
| `writeoutangle` | Write angle terms from database | `True` |
| `writeoutpolarize` | Write polarise terms from database | `True` |
| `writeouttorsion` | Write torsion terms from database | `True` |
| `inputkeyfile` | Append parameters from this key file | — |
| `quickdatabasesearch` | Skip exhaustive SMARTS re-matching | `False` |

---

### Torsion Parameterisation

| Keyword | Description | Default |
|---------|-------------|---------|
| `rotalltors` | Scan all rotatable bonds (ignore database matches) | `False` |
| `toroptmethod` | QM method for torsion geometry optimisation | `wB97X-D` |
| `toroptbasisset` | Basis set for torsion geometry optimisation | `6-31G*` |
| `toroptbasissetfile` | Basis set file for torsion optimisation | `6-31g_st_.0.gbs` |
| `torspbasisset` | Basis set for torsion single-point | `6-311+G*` |
| `torspbasissetfile` | Basis set file for torsion single-point | `6-311+g_st_.0.gbs` |
| `torspbasissethalogen` | Torsion single-point basis set for halogenated molecules | — |
| `toroptmethodlist` | Comma-separated per-job optimisation method overrides | — |
| `torspmethodlist` | Comma-separated per-job single-point method overrides | — |
| `tortor` | Enable torsion-torsion (2D) coupling | `False` |
| `dontrotbndslist` | Comma-separated bond pairs to exclude from scanning | — |
| `onlyrottortorlist` | Comma-separated bond pairs restricted to tor-tor scanning | — |
| `nonaroringtor1dscan` | Enable 1D scanning of non-aromatic ring torsions | `False` |
| `xtbtorresconstant` | XTB torsion restraint force constant | 5.0 |
| `maxtorresnitrogen` | Maximum torsion grid points for N-containing bonds | 2 |
| `defaultmaxtorsiongridpoints` | Maximum torsion grid points per bond | 40 |

---

### Fragmentation

| Keyword | Description | Default |
|---------|-------------|---------|
| `smallmoleculefragmenter` | Use small-molecule (non-WBO) fragmentation | `False` |
| `quickdatabasesearch` | Skip exhaustive re-matching in fragment database | `False` |

---

### Van der Waals Parameterisation

| Keyword | Description | Default |
|---------|-------------|---------|
| `dovdwscan` | Refine vdW parameters with ab initio dimer data | `False` |
| `onlyvdwatomlist` | Comma-separated atom indices to include in refinement | — |
| `fixvdwtyperadii` | Comma-separated vdW types whose radii are fixed | — |
| `accuratevdwsp` | Use CBS extrapolation (aug-cc-pV[TQ]Z, Psi4 only) | `False` |
| `vdwprmtypestofit` | vdW parameter types to fit (`S`=sigma, `T`=taper) | `S,T` |
| `homodimers` | Use homo-molecular dimers for vdW fitting | `False` |
| `vdwmaxqmstartingpointspertype` | Max QM starting points per vdW type | 1 |
| `vdwmaxtinkergridpoints` | Max Tinker grid points per vdW probe | 50 |

---

### Iodine / Heavy-Halogen Overrides

When iodine or other heavy halogens are present, the default small basis sets
are replaced automatically.  These keywords allow manual control:

| Keyword | Description | Default |
|---------|-------------|---------|
| `iodineoptbasisset` | Basis set for opt of iodine-containing molecules | `def2-svp` |
| `iodineoptbasissetfile` | Basis set file for iodine optimisation | `def2-svp.1.gbs` |
| `iodinedmabasissetfile` | Basis set file for iodine DMA | `def2-svp.1.gbs` |
| `iodineespbasissetfile` | Basis set file for iodine ESP | `def2-tzvpp.1.gbs` |
| `iodinetoroptbasisset` | Basis set for iodine torsion optimisation | `def2-svp` |
| `iodinetoroptbasissetfile` | Basis set file for iodine torsion opt | `def2-svp.1.gbs` |
| `iodinetorspbasisset` | Basis set for iodine torsion single-point | `def2-svp` |
| `iodinetorspbasissetfile` | Basis set file for iodine torsion SP | `def2-svp.1.gbs` |

---

### ML Potentials

| Keyword | Description | Default |
|---------|-------------|---------|
| `anienvname` | Conda environment name for ANI potentials | `ani` |
| `xtbenvname` | Conda environment name for XTB | `xtbenv` |
| `anifmax` | ANI force convergence threshold | 0.05 |
| `xtbmethod` | XTB method level (1 or 2) | 2 |
| `fennixmodelname` | FeNNol model name | `ani2x_model0` |

---

### Pipeline Control

| Keyword | Description | Default |
|---------|-------------|---------|
| `dryrun` | Validate inputs without running QM | `False` |
| `databasematchonly` | Stop after database matching | `False` |
| `setupfragjobsonly` | Stop after setting up fragment jobs | `False` |
| `optonly` | Stop after geometry optimisation | `False` |
| `checkinputonly` | Validate inputs and exit | `False` |
| `generateinputfilesonly` | Generate QM input files and exit | `False` |
| `genprotstatesonly` | Generate protonation states and exit | `False` |

> **Note:** `--dry-run` on the CLI sets `dryrun=True` automatically.
