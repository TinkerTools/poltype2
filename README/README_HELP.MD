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

* All booleans accept `True` or `False` values.
* Defaults are chosen to maximise accuracy and efficiency.
* Large molecules may require extra memory for QM jobs.

---

### CLI Options

```
python -m poltype [OPTIONS]

  -c, --config FILE        poltype.ini path (default: poltype.ini)
  -w, --work-dir DIR       Working directory (default: current directory)
  -v, --verbose            Enable debug logging
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
| `numproc` | Number of processors | 80 % of available |
| `maxmem` | Maximum RAM | 80 % of available |
| `maxdisk` | Maximum scratch disk | 80 % of available |
| `consumptionratio` | Fraction of node resources to use | 0.8 |
| `coresperjob` | Cores per parallel QM job | 2 |
| `jobsatsametime` | Maximum concurrent fragment jobs | `numproc / coresperjob` |
| `parentjobsatsametime` | Number of Poltype jobs on same host | 1 |

---

### Protonation State

| Keyword | Description | Default |
|---------|-------------|---------|
| `addhydrogens` | Add hydrogens to unprotonated atoms | `False` |
| `addhydrogentocharged` | Add hydrogens to under-valent C/N | `True` |
| `genprotstatesonly` | Output protonation states and exit | `False` |

---

### Quantum Chemistry Packages

| Keyword | Description | Default |
|---------|-------------|---------|
| `use_gaus` | Use Gaussian for all QM calculations | `False` |
| `use_gausoptonly` | Use Gaussian only for geometry optimisation | `False` |

By default, PySCF is used.  Set `use_gaus=True` for Gaussian or leave
default for Psi4/PySCF.

---

### Atom Typing

| Keyword | Description | Default |
|---------|-------------|---------|
| `atmidx` | Starting type number for atom types | 401 |
| `usesymtypes` | Use symmetry to reduce and average atom types | `True` |
| `indextotypefile` | Custom atom-index-to-type mapping file | `None` |

---

### Geometry Optimisation

| Keyword | Description | Default |
|---------|-------------|---------|
| `optbasisset` | Basis set for geometry optimisation | `6-31G*` |
| `optmethod` | QM method for geometry optimisation | `MP2` |
| `optpcm` | Use PCM for optimisation (`True`/`False`/`auto`) | `auto` |
| `pcm_auto` | Auto-activate PCM for zwitterions | `True` |
| `generateextendedconf` | Use extended conformation with restrained torsions | `True` |
| `userconformation` | Use user-supplied conformation | `False` |
| `optloose` | Use `GAU_LOOSE` convergence criteria | `False` |
| `optconvergence` | Convergence level (`LOOSE`/`VERYLOOSE`) | `LOOSE` |

---

### Electrostatics (DMA and ESP)

| Keyword | Description | Default |
|---------|-------------|---------|
| `dmamethod` | QM method for DMA | `MP2` |
| `dmabasisset` | Basis set for DMA | `6-311G**` |
| `espmethod` | QM method for ESP fitting | `MP2` |
| `espbasisset` | Basis set for ESP fitting | `aug-cc-pvtz` |
| `sameleveldmaesp` | Share a single QM job for DMA and ESP | `False` |
| `numespconfs` | Number of conformations for ESP fitting | 1 |
| `espextraconflist` | Extra conformer files for ESP fitting | `[]` |
| `usepoleditframes` | Use POLEDIT frames instead of Poltype detection | `True` |
| `new_gdma` | Use grid-based DMA algorithm (experimental) | `False` |

---

### Database Matching

| Keyword | Description | Default |
|---------|-------------|---------|
| `writeoutmultipole` | Write multipole terms from database | `True` |
| `writeoutbond` | Write bond terms from database | `True` |
| `writeoutangle` | Write angle terms from database | `True` |
| `writeoutpolarize` | Write polarise terms from database | `True` |
| `writeouttorsion` | Write torsion terms from database | `True` |
| `inputkeyfile` | Append parameters from this key file | `None` |
| `prmmodfile` | Comma-separated parameter-modification patch files | `None` |

---

### Torsion Parameterisation

| Keyword | Description | Default |
|---------|-------------|---------|
| `dontdotor` | Skip torsion parameterisation | `False` |
| `rotalltors` | Scan all rotatable bonds (ignore database) | `False` |
| `foldnum` | Fourier terms in torsion fit | 6 |
| `toroptmethod` | Optimisation method for torsion scanning (`xtb`/`ANI`/`AMOEBA`) | `xtb` |
| `torspmethod` | Single-point method for torsion energies | `wB97XD` |
| `toroptbasisset` | Basis set for torsion optimisation | `6-31G*` |
| `torspbasisset` | Basis set for torsion single points | `6-311+G*` |
| `torsppcm` | PCM for torsion SP (`True`/`False`/`auto`) | `auto` |
| `toroptpcm` | PCM for torsion optimisation (`True`/`False`/`auto`) | `auto` |
| `torsionrestraint` | Restraint force constant (kcal/mol/rad^2) | 1641.4 |
| `tortor` | Enable torsion-torsion coupling | `False` |
| `onlyrotbndslist` | Only rotate these bonds (e.g. `1 2, 3 4`) | — |
| `dontrotbndslist` | Skip rotation of these bonds | — |

---

### Fragmentation

| Keyword | Description | Default |
|---------|-------------|---------|
| `dontfrag` | Disable fragmentation | `False` |
| `WBOtol` | Wiberg bond-order tolerance for parent/fragment match | 0.05 |
| `maxgrowthcycles` | Maximum fragment growth iterations | 4 |

---

### Van der Waals Parameterisation

| Keyword | Description | Default |
|---------|-------------|---------|
| `dovdwscan` | Refine vdW parameters with ab initio dimer data | `False` |
| `onlyvdwatomlist` | vdW atom indices to include in refinement | — |
| `fixvdwtyperadii` | vdW types whose radii are fixed during fitting | — |
| `accuratevdwsp` | Use CBS extrapolation with aug-cc-pV[TQ]Z | `False` |

---

### Pipeline Control

| Keyword | Description | Default |
|---------|-------------|---------|
| `dry_run` | Validate inputs without running QM | `False` |
| `skip_validation` | Skip the validation stage | `False` |
| `database_match_only` | Stop after database matching | `False` |
| `setup_frag_jobs_only` | Stop after setting up fragment jobs | `False` |
| `opt_only` | Stop after geometry optimisation | `False` |
| `check_input_only` | Validate inputs and exit | `False` |
