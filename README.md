# Poltype2: Automation of AMOEBA Polarizable Force Field for Small Molecules

## Objective

Given an input chemical structure, all force-field parameters are automatically
assigned from a database or derived via fitting to ab initio data generated on
the fly.

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![License][license-shield]][license-url]

[contributors-shield]: https://img.shields.io/github/contributors/TinkerTools/poltype2.svg?style=for-the-badge
[contributors-url]: https://github.com/TinkerTools/poltype2/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/TinkerTools/poltype2.svg?style=for-the-badge
[forks-url]: https://github.com/TinkerTools/poltype2/network/members
[stars-shield]: https://img.shields.io/github/stars/TinkerTools/poltype2.svg?style=for-the-badge
[stars-url]: https://github.com/TinkerTools/poltype2/stargazers
[issues-shield]: https://img.shields.io/github/issues/TinkerTools/poltype2.svg?style=for-the-badge
[issues-url]: https://github.com/TinkerTools/poltype2/issues
[license-shield]: https://img.shields.io/github/license/TinkerTools/poltype2.svg?style=for-the-badge
[license-url]: https://github.com/TinkerTools/tinker/blob/release/LICENSE.pdf

## Please Cite

Walker, B., Liu, C., Wait, E., Ren, P., J. Comput. Chem. 2022, 1.
<https://doi.org/10.1002/jcc.26954>

Wu JC, Chattree G, Ren P. Automation of AMOEBA polarizable force field
parameterization for small molecules. Theor Chem Acc. 2012;131(3):1138.
doi:10.1007/s00214-012-1138-6

## License

[Tinker License](https://github.com/TinkerTools/tinker/blob/release/LICENSE.pdf)

---

## Features

### Input Features

* Automated total charge assignment from bond orders and element valences
* Dominant ionization state enumeration (pH 7) and tautomer enumeration
* Smart memory / CPU resource defaults for QM jobs

### Parameterization Features

* 10-stage modular pipeline with checkpoint/resume support
* Molecule fragmentation to speed up QM calculations
* Parallelised QM job submission (Psi4, Gaussian, PySCF)
* SMARTS-based atom typing and database parameter matching
* QM dimer data generation for van der Waals fitting
* Torsion-Torsion coupling
* Expanded torsion parameter database
* Post-run validation (energy convergence, QM/MM RMSD)

---

## Quick Start

### Installation

```bash
pip install -e ".[dev]"
```

See [Installation Guide](README/README_INSTALL.MD) for full details on
installing Tinker, GDMA, and quantum chemistry packages.

### Minimum Usage

Create a `poltype.ini` file:

```ini
structure=methylamine.sdf
```

Run Poltype:

```bash
python -m poltype --config poltype.ini --work-dir ./output
```

The resulting `final.xyz` and `final.key` are the structure and parameter
files you need.

### CLI Options

```
python -m poltype --help

Options:
  -c, --config       Path to poltype.ini (default: poltype.ini)
  -w, --work-dir     Working directory for output (default: current dir)
  -v, --verbose      Enable verbose/debug logging
  -q, --quiet        Suppress all output except warnings and errors
  --resume           Resume from the last checkpoint
  --no-checkpoint    Disable checkpoint saving
  --dry-run          Validate inputs without running QM jobs
  --version          Show version and exit
```

---

## Pipeline Architecture

Poltype2 uses a modular 10-stage pipeline:

| # | Stage | Description |
|---|-------|-------------|
| 1 | **Input Preparation** | Load molecule, assign charges, validate structure |
| 2 | **Geometry Optimisation** | QM geometry optimisation (MP2/6-31G*) |
| 3 | **ESP Fitting** | Electrostatic potential grid computation |
| 4 | **Multipole** | Distributed multipole analysis (GDMA) |
| 5 | **Atom Typing** | SMARTS-based atom type assignment |
| 6 | **Database Match** | Look up existing parameters |
| 7 | **Fragmentation** | Fragment large molecules for torsion scanning |
| 8 | **Torsion Fitting** | Torsion angle scanning and parameter fitting |
| 9 | **Validation** | Energy convergence and QM/MM RMSD checks |
| 10 | **Finalisation** | Write `final.xyz`, `final.key`, and summary report |

Each stage is independently testable and communicates via an artifact
dictionary in the pipeline context.

---

## Repository Layout

```
poltype2/
├── poltype/                 # Main Python package
│   ├── __main__.py          # CLI entry point
│   ├── config/              # Configuration schema and INI loader
│   ├── molecule/            # RDKit-based molecule representation
│   ├── qm/                  # QM backend abstraction (Psi4, Gaussian, PySCF)
│   ├── typing/              # SMARTS-based atom typing
│   ├── database/            # Parameter database lookup
│   ├── fragmentation/       # Rule-based and WBO molecular fragmentation
│   ├── validation/          # Post-run energy and RMSD validation
│   ├── output/              # Output writers (XYZ, KEY, summary)
│   ├── pipeline/            # Pipeline framework (stages, runner, events, checkpoints)
│   └── errors.py            # Custom exception hierarchy
├── ParameterFiles/          # AMOEBA parameter databases and SMARTS mappings
├── BasisSets/               # Quantum chemistry basis set files
├── FennixModels/            # ANI machine-learning potentials
├── Examples/                # Example parameterisation runs
├── tests/unit/              # Unit tests (597 tests)
├── README/                  # Extended documentation
└── pyproject.toml           # Python packaging metadata
```

---

## Documentation

| Document | Description |
|----------|-------------|
| [Installation Guide](README/README_INSTALL.MD) | Tinker, GDMA, QM packages, and Python environment setup |
| [Advanced Usage](README/README_HELP.MD) | All `poltype.ini` keywords and CLI options |
| [Output Files](README/README_OUTPUT.MD) | Description of generated files (`final.xyz`, `final.key`, logs) |
| [Protocol Details](README/README_PROTOCOL.MD) | Technical details of the parameterisation workflow |
| [File Manifest](README/README_MANIFEST.MD) | Intermediate and output file descriptions |

---

## Input Preparation

* Input structure should be a file with bond order information (`.sdf`,
  `.mol`, `.mol2`, etc.)
* Cartesian XYZ, Tinker XYZ, and PDB are not recommended—these will be
  converted to SDF and bond orders guessed from atomic distances
* 2D structures are automatically converted to 3D coordinates
* Total charge is computed from formal charges on each atom
* Use `addhydrogens=True` to add missing hydrogens
* Use `addhydrogentocharged=False` to prevent protonation of charged atoms

## Validation Checks

After parameterisation, the `ValidationStage` runs `EnergyValidator` checks:

* Geometry optimisation energy is finite and the run converged
* Torsion energy range is within 50 kcal/mol (configurable)
* QM vs MM torsion RMSD is below 3 kcal/mol (configurable)

Validation is skipped automatically when `config.skip_validation = True`
or when `--dry-run` is active.

