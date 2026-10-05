# Poltype2: Automation of AMOEBA Polarizable Force Field for Small Molecules

## Objective
Given an input chemical structure, all parameters can be automatically assigned from a database or derived via fitting to ab initio data generated on the fly.

[![forthebadge made-with-python](http://ForTheBadge.com/images/badges/made-with-python.svg)](https://www.python.org/)

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]


[contributors-shield]: https://img.shields.io/github/contributors/TinkerTools/poltype2.svg?style=for-the-badge
[contributors-url]: https://github.com/TinkerTools/poltype2/forks/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/TinkerTools/poltype2.svg?style=for-the-badge
[forks-url]: https://github.com/TinkerTools/poltype2/network/members
[stars-shield]: https://img.shields.io/github/stars/TinkerTools/poltype2.svg?style=for-the-badge
[stars-url]: https://github.com/TinkerTools/poltype2/stargazers
[issues-shield]: https://img.shields.io/github/issues/TinkerTools/poltype2.svg?style=for-the-badge
[issues-url]: https://github.com/TinkerTools/poltype2/issues
[license-shield]: https://img.shields.io/github/license/TinkerTools/poltype2.svg?style=for-the-badge
[license-url]: https://github.com/TinkerTools/tinker/blob/release/LICENSE.pdf


## Please Cite

Walker, B., Liu, C., Wait, E., Ren, P., J. Comput. Chem. 2022, 1. https://doi.org/10.1002/jcc.26954

Wu JC, Chattree G, Ren P. Automation of AMOEBA polarizable force field parameterization for small molecules. Theor Chem Acc. 2012 Feb 26;131(3):1138. doi: 10.1007/s00214-012-1138-6. PMID: 22505837; PMCID: PMC3322661.

## License
[Tinker License](https://github.com/TinkerTools/tinker/blob/release/LICENSE.pdf)

## 📚 Documentation Overview

* Please read 👇🙏

### Features

* Parameterization Input Features
    * Automated total charge assignment
    * Dominant ionization state enumeration (pH=7) / tautomer enumeration
    * Smart memory resource defaults for QM jobs
* Parameterization Features
    * Molecule fragmenter to speed up QM calculations
    * Parallelized job submission for QM jobs
    * Psi4/Gaussian quantum package support
    * QM dimer data generation for vdW fitting
    * Torsion-Torsion coupling
    * Expanded torsion database


[💻 Program Installation](README/README_INSTALL.MD)


### Automated AMOEBA Ligand Parameterization (Main function of Poltype)

* [Parameterization Input Preparation](#parameterization-input-preparation)
* [Ligand Protonation State Generation](#ligand-protonation-state-generation)
* [Minimum Example Usage Parameterization](#minimum-example-usage-parameterization)
* [Default Resource Consumption](#default-resource-consumption)
* [Parameterization Sanity Checks](#parameterization-sanity-checks)
* [Keyword Reference](#keyword-reference) — every `poltype.ini` keyword recognized by `poltype.py`
* [💻 Advanced Program Usage](README/README_HELP.MD) — curated walkthroughs and recommended-input recipes
* [Parameterization Output Files](README/README_OUTPUT.MD)
    * [Final XYZ File](README/README_OUTPUT.MD)
    * [Final Key File](README/README_OUTPUT.MD)
    * [Poltype Log File](README/README_OUTPUT.MD)
    * [OPENME Plots](README/README_OUTPUT.MD)
* [Parameterization Examples](Examples/Parameterization)
* [Automated AMOEBA Ligand Parameterization How It Works](README/README_PROTOCOL.MD)


---------------------------------------------------------------------------------------------


### Parameterization Input Preparation
* The input structure can be given as a filetype with bond order information (sdf, mol, mol2, etc.)
* Cartesian XYZ, tinker XYZ and PDB are not recommended; they will be converted to an SDF file and bond orders will be guessed based on atomic distances.
* If a 2D structure is given, Poltype will generate 3D coordinates for you.
* Total charge is determined by computing the formal charge for each atom.
* Formal atom charge will be assigned via the input number of bonds and bond order for surrounding bonds and element of each atom.
* Optional keywords exist to add missing hydrogens.
* If carbon or nitrogen are negatively charged, then Poltype will assume you meant to have hydrogens on those atoms and protonate to neutral charge. For all other elements, if there exists a formal charge (due to input bond order and valence electrons of element), then Poltype will detect the formal charge for you and compute total charge from sum of each formal atom charge. A warning message is printed when hydrogen atoms are added.
* Special radical charge states require additional information in the input file specifying which atom is a radical.
* If you do not want charged atoms to be protonated then use ``addhydrogentocharged=False``.

### Ligand Protonation State Generation
* Dominant ionization states at pH 7 are enumerated and SDF files are generated via Dimorphite-DL (`IonizationState_0.sdf`, `IonizationState_1.sdf`, ...).
* Tautomer states are enumerated via rdkit and the first tautomer is the canonical tautomer `TautomerState_0.sdf`.
* Use ``genprotstatesonly`` to quit the program after generating dominant ionization states at pH=7 and tautomers.

### Minimum Example Usage Parameterization

__All input arguments are specified in a `poltype.ini` file__
```
structure=methylamine.sdf
```
* Navigate to the directory containing `poltype.ini` and the `.sdf` file, and run:

```shell
nohup python /path_to_poltype/poltype.py &
```

`final.xyz` and `final.key` are the resulting structure and parameter files you will need.
* After Poltype finishes, check the ``OPENME`` folder for torsion-fitting and ESP-fitting results.

### Default Resource Consumption
* By default, Poltype computes the number of fragment Poltype jobs (or any QM job if the fragmenter is not being used) to run in parallel as the floor of the input number of cores divided by the number of cores per job (default of 2).
* RAM, cores, and disk space can all be detected and a consumption ratio of 80% is used by default. Fragment RAM, disk, and cores are divided evenly by the number of Poltype jobs in parallel.

### Parameterization Sanity Checks
* MM = Molecular Mechanics (AMOEBA model), QM = Quantum Mechanics
* Check for 2D coordinates and generate 3D coordinates at the beginning of the program.
* Check to ensure refined multipoles enable MM to model the QM potential grid well; raise an error if not.
* Check to ensure QM and MM dipoles are very similar; raise an error if not.
* Check to ensure the final minimized MM structure is similar to the QM-optimized geometry; raise an error if not.
* Check for any missing van der Waals parameters at the end of the program; raise an error.
* Check for any missing multipole parameters at the end of the program; raise an error.
* Check for any zeroed-out torsion parameters in the final key file; check if the fragmenter is transferring torsions properly.


---------------------------------------------------------------------------------------------


## Keyword Reference

This section describes every keyword recognized by `poltype.py` when reading a `poltype.ini` configuration file. Each entry lists the keyword's purpose, accepted value type, and default.

**Conventions**
- Booleans accept `True` / `False` (case-insensitive). Writing the keyword alone (no `=`) is treated as `True`.
- Three-state PCM switches (`optpcm`, `toroptpcm`, `torsppcm`) additionally accept `auto` (encoded internally as `-1`).
- List keywords are comma-separated; index lists use space-separated entries inside each comma-separated element (e.g. `onlyrotbndslist=1 2, 3 4`).
- "Vestigial" means the keyword is parsed by `poltype.py` but no downstream code reads the resulting attribute; setting it has no effect.

**Keyword categories**
1. [Required Input](#1-required-input)
2. [Resource and Job Control](#2-resource-and-job-control)
3. [Run Mode / Workflow Gates](#3-run-mode--workflow-gates)
4. [Protonation State](#4-protonation-state)
5. [QM Package Selection](#5-qm-package-selection)
6. [Atom Typing](#6-atom-typing)
7. [Geometry Optimization](#7-geometry-optimization)
8. [Distributed Multipole Analysis (DMA) and ESP](#8-distributed-multipole-analysis-dma-and-esp)
9. [Database Matching / Pre-fit Key Output](#9-database-matching--pre-fit-key-output)
10. [Torsion Parameterization](#10-torsion-parameterization)
11. [Fragmenter](#11-fragmenter)
12. [Van der Waals Parameterization](#12-van-der-waals-parameterization)
13. [QM Output File Names](#13-qm-output-file-names)
14. [Tolerances and Error Reporting](#14-tolerances-and-error-reporting)
15. [Legacy / External-Workflow Keywords](#15-legacy--external-workflow-keywords)


### 1. Required Input

#### `structure`
SDF input molecule for parameterization. Must contain bond orders consistent with the desired total charge. **Required.** Default: `None`.

#### `totalcharge`
Integer net charge of the molecule. If unset, Poltype infers it from the input SDF. Default: `None`.


### 2. Resource and Job Control

#### `numproc`
Number of processor cores Poltype is allowed to use. Default: ~80% of node cores.

#### `maxmem`
Maximum RAM available to QM jobs (e.g. `20GB`). Default: ~80% of node memory.

#### `maxdisk`
Maximum scratch disk allotted to QM jobs (e.g. `100GB`). Default: ~80% of free scratch.

#### `consumptionratio`
Float in `(0,1]` controlling what fraction of node resources Poltype auto-claims when `numproc` / `maxmem` / `maxdisk` are unset. Default: `0.8`.

#### `coresperjob`
Number of cores per parallel sub-job (e.g. one fragment QM job). Default: `2`.

#### `parentjobsatsametime`
Number of top-level Poltype jobs sharing the same host. Poltype divides resources by this number. Default: `1`.

#### `lastlogfileupdatetime`
Hours a log file must be idle before Poltype assumes the associated job has died. Default: `1`.

#### `sleeptime`
Polling interval (seconds) between job-status checks. Default: `0.1`.

#### `bashrcpath`
Absolute path to a bash startup script Poltype sources before invoking external QM engines. Default: `None`.


### 3. Run Mode / Workflow Gates

#### `checkinputonly`
If `True`, validate the input and exit without running any QM. Default: `False`.

#### `generateinputfilesonly`
If `True`, write all QM input files but do not submit jobs. Default: `False`.

#### `setupfragjobsonly`
If `True`, only prepare fragment subdirectories/inputs, then stop. Default: `False`.

#### `databasematchonly`
If `True`, perform only database parameter assignment and exit (no QM, no fitting). Default: `False`.

#### `qmonly`
If `True`, run all QM calculations and stop before final parameter fitting. Default: `False`.

#### `firstoptfinished`
Tells Poltype the initial geometry optimization is already complete; resume from DMA. Default: `False`.

#### `isfragjob`
Marker that the current run is a fragment sub-job (set internally by Poltype). Default: `False`.

#### `fragmenterdebugmode`
Enable verbose debugging in the fragmenter. Default: `False`.

#### `toroptdebugmode`
Enable verbose debugging during torsion optimization. Default: `False`.

#### `tordebugmode`
Enable verbose debugging during torsion fitting. Default: `False`.

#### `debugmode`
Enable verbose top-level debugging in the parameterization driver. Matched only when none of `fragmenterdebugmode` / `tordebugmode` / `toroptdebugmode` is on the line. Default: `False`.

#### `printoutput`
Stream QM/Tinker tool output to stdout in addition to log files. Default: `False`.

#### `deleteallnonqmfiles`
Clean up auxiliary scratch (non-QM output) when finished. Default: `False`.

#### `jobsatsametime`
Maximum number of QM sub-jobs running in parallel on the host. `0` (default) means auto-compute as `floor(numproc/coresperjob)`, capped at `maxjobsatsametime=10`. Default: `0`.

#### `optonly`
Run only the initial geometry optimization and stop. The match excludes lines containing `gaus` so it does not collide with `use_gausoptonly`. Default: `False`.

#### `useuniquefilenames`
Vestigial — assigned in the ini reader and never read elsewhere. Default (if set bare): `True`.

#### `checktraj`
Vestigial — assigned in the ini reader and never read elsewhere. Default (if set bare): `True`.


### 4. Protonation State

#### `addhydrogens`
Add hydrogens to unprotonated heavy atoms before parameterization. Default: `False`.

#### `addhydrogentocharged`
Add hydrogens to charged-looking atoms (C with valence < 4, N with valence 2). Default: `True`.

#### `genprotstatesonly`
Write protonation states at pH 7 to `ProtonationState_*.mol` files and exit. Default: `False`.


### 5. QM Package Selection

#### `use_gaus`
Use Gaussian for **all** QM steps (opt, DMA, ESP, torsion). Matched only when `opt` is not in the line so it does not collide with `use_gausoptonly`. Default: `False`.

#### `use_gausoptonly`
Use Gaussian only for the initial geometry optimization; use Psi4 for everything else. Default: `False`.

#### `use_gauPCM`
Use Gaussian's PCM implementation rather than Psi4's. Default: `False`.

#### `use_gau_vdw`
Use Gaussian for vdW dimer single-point energies. Default: `False`.

#### `use_psi4_geometric_opt`
Use Psi4's GeomeTRIC engine for optimization. Default: `True`.

#### `dont_use_pyscf`
Disable the PySCF backend. Default: `False`.

#### `use_qmopt_vdw`
Use the post-QM-opt geometry (not the AMOEBA-minimized one) when sampling vdW dimers. Default: `False`.

#### `scfmaxiter`
Maximum SCF iterations for QM single points. Default: `500`.

#### `optmaxcycle`
Maximum geometry-optimization steps. Default: `400`.

#### `gausoptcoords`
Gaussian-style coordinate keyword inserted into the optimization route line (e.g. `cartesian`). Default: empty.

#### `pyscf_opt_meth`
Exchange-correlation functional used by the PySCF optimizer (written into the generated PySCF input as the `xc` argument). Default: `wb97x-d3bj`.

The default requires PySCF's dispersion support (`pyscf-dispersion`; see [README_INSTALL](README/README_INSTALL.MD)). The older default `wb97x_d3` only resolves on **pyscf <= 2.7.0**, where the underscore spelling slipped past the dispersion blacklist and mapped to libxc `HYB_GGA_XC_WB97X_D3` (XC part only, no additive D3 term). pyscf 2.8+ rejects both spellings, so use `wb97x-d3bj` (wB97X-V XC + D3BJ) there. The two are different functionals — energies and optimized geometries will not match pre-2.8 `wb97x_d3` results.

#### `pyscf_sol_imp`
PySCF implicit-solvent model: `IEF-PCM`, `C-PCM`, `SS(V)PE`, or `COSMO`. Default: `IEF-PCM`.

#### `pyscf_sol_eps`
PySCF solvent dielectric constant. Default: `78.3553` (water).


### 6. Atom Typing

#### `atmidx`
Alias of `prmstartidx` — starting number for assigned AMOEBA type indices. Default: `401`.

#### `prmstartidx`
Starting AMOEBA type index. Default: `401`.

#### `usesymtypes`
Use molecular symmetry to merge equivalent atom types. Default: `True`.

#### `indextotypefile`
Path to a custom file mapping 1-indexed atom numbers to target type numbers (one mapping per line). Default: `None`.

#### `indextompoleframefile`
File providing multipole-frame definitions (in the same format as the first key file) instead of Poltype's automatic detection. Default: `None`.

#### `forcefield`
Force-field label written into the resulting key file. Default: `AMOEBA`.

#### `paramhead`
Path to the header parameter file appended to the top of the output key file. Default: `ParameterFiles/amoebabio18_header.prm`.

#### `prmfilepath`
Absolute path to the base AMOEBA `.prm` library used for typing. Default: `ParameterFiles/amoebabio18.prm`.

#### `parentname`
Override the parent-molecule label used when this run is a fragment of a larger molecule. Default: derived from input filename.

#### `poltypepath`
Override the Poltype installation path used internally. Matched only when `poltypepathlist` is not on the line. Default: directory containing `poltype.py`.

#### `polareps`
Vestigial — assigned in the ini reader and never read elsewhere. Default: unset.


### 7. Geometry Optimization

#### `optmethod`
QM method for geometry optimization. Default: `MP2`.

#### `optbasisset`
Basis set for geometry optimization. Default: `6-31G*`.

#### `optpcm`
Use PCM for optimization. `True` / `False` / `auto`. `auto` enables PCM for zwitterions. Default: `auto` (encoded as `-1`).

#### `pcm_auto`
Master switch for automatic PCM activation on zwitterions. Default: `True`.

#### `dontusepcm`
Logical inverse of `pcm_auto` — disables automatic PCM. Default: `False`.

#### `optloose`
Use Psi4 `GAU_LOOSE` convergence criteria for opt and torsion opt. Default: `True`.

#### `optconvergence`
Convergence preset (uppercased on read). Values handled in `optimization.py` / `torsiongenerator.py`: `LOOSE` (GAU_LOOSE criteria), `VERYLOOSE` (loose gradient only), `NULL` (default Psi4/Gaussian criteria). Default: `LOOSE`.

#### `freq`
Run a vibrational-frequency calculation after optimization. Default: `False`.

#### `generateextendedconf`
Generate an extended starting conformer (torsions trans-restrained) for opt. Default: `True`.

#### `userconformation`
Use the conformer from the input SDF directly with torsions restrained (skips extended-conf generation). Default: `False`.

#### `userxyzgeometry`
Path to a Tinker XYZ already-optimized geometry; Poltype skips opt and proceeds to DMA. Default: empty.

#### `generate_symm_frag_conf`
Generate symmetric-fragment conformers during conformer expansion. Default: `False`.

#### `sp2aniline`
Force planar `c-NH2` (sp2 anilines); also assigns `anglep` / `opbend`. Used with `use_gaus`. Default: `False`.

#### `nonplanarphenol`
Force non-planar phenol geometry. Default: `False`.

#### `force_anti_acid`
Put carboxylic acid `O=C-O-H` in the anti (180°) conformation when generating the extended conformer; the QM optimization keeps it frozen, so multipoles are derived from the anti conformer. Requires `generateextendedconf=True`. Default: `False`.


### 8. Distributed Multipole Analysis (DMA) and ESP

#### `dmamethod`
QM method for the DMA single point. Default: `MP2`.

#### `dmabasisset`
Basis set for the DMA single point. Default: `6-311G**`.

#### `espmethod`
QM method for ESP fitting. Default: `MP2`.

#### `espbasisset`
Basis set for ESP fitting. Default: `aug-cc-pVTZ`.

#### `popbasisset`
Basis set used for the Gaussian population-analysis (`Pop=SaveMixed`) step in ESP setup. No class-level default — must be set in `poltype.ini` if that code path is used.

#### `espfit`
Fit electrostatic parameters to the QM ESP grid. Default: `True`.

#### `esprestweight`
Restraint weight (written as `RESP-WEIGHT`) applied to derived multipoles relative to the DMA initial guess during ESP fitting. Default: `1`.

#### `espgrad`
Convergence-gradient setting passed to the Tinker `potential` executable during ESP fitting. Default: `0.1`.

#### `espextraconflist`
Comma-separated list of extra conformer files (SDF or equivalent, same atom order) to add to the ESP fit. Default: empty.

#### `numespconfs`
Vestigial — assigned in the ini reader and never read elsewhere. The multi-conformer ESP-fitting path appears to have been replaced by `espextraconflist`. Default: `1`.

#### `sameleveldmaesp`
Run QM only once for DMA and reuse the wavefunction for ESP. Default: `False`.

#### `adaptiveespbasisset`
Downgrade `aug-cc-pVTZ` to `aug-cc-pVDZ` when the molecule has ≥20 heavy atoms. Default: `False`.

#### `new_gdma`
Use the new grid-based DMA4 algorithm (experimental). Also auto-enables `scaleandfixdipole`. Default: `False`.

#### `gdmacommand_*`
Free-form prefix to inject additional GDMA commands. `gdmacommand_Radius_S=0.80` becomes `Radius S 0.80`. Default: none.

#### `scaleandfixdipole`
For DMA4: scale and fix dipole components flagged as overly large (e.g. COO atoms). Default: `False` (auto-enabled by `new_gdma=True`).

#### `scalebigmultipole`
When DMA produces unusually large multipoles, scale them down before high-level ESP fitting. Default: `False`.

#### `fragbigmultipole`
When DMA produces unusually large multipoles, derive replacements from a smaller fragment. Default: `True`.

#### `atomidsfordmafrag`
Comma-separated atom indices to treat as "Big-ID multipole" atoms with `fragbigmultipole`. Default: empty.

#### `usepoleditframes`
Use POLEDIT-style frame definitions instead of Poltype's frame detection. Default: `True`.

#### `potentialoffset`
Offset (written as `potential-offset`) applied when matching grid potentials. Default: `1.0`.

#### `chargethreshold`
Absolute magnitude above which an atomic monopole is flagged as a "big multipole" candidate for scaling or fragment-based replacement. Only C and N atoms are considered. Default: `1.5`.

#### `dipolethreshold`
Absolute magnitude above which any atomic dipole component is flagged as a "big multipole". Default: `1.5`.

#### `quadrupolethreshold`
Absolute magnitude above which any atomic quadrupole component is flagged as a "big multipole". Default: `2.4`.

#### `fitqmdipole`
Add the QM molecular dipole as an additional restraint during ESP fitting. Default: `False`.

#### `dipoletol`
Relative tolerance (`|qm-mm|/qm`) for the molecular-dipole comparison. Default: `0.5` (i.e. 50%).

#### `absdipoletol`
Absolute tolerance (D) on `|qm-mm|`. An error is raised only when **both** the relative (`dipoletol`) and absolute (`absdipoletol`) gates are exceeded. Default: `0.5`.

#### `suppressdipoleerr`
Demote dipole-mismatch errors to warnings. Default: `False`.


### 9. Database Matching / Pre-fit Key Output

#### `writeoutmultipole`
Write database-matched multipole terms into the pre-fit key file. Default: `True`.

#### `writeoutbond`
Write database-matched bond terms. Default: `True`.

#### `writeoutangle`
Write database-matched angle terms. Default: `True`.

#### `writeoutpolarize`
Write database-matched polarize terms. Default: `True`.

#### `writeouttorsion`
Write database-matched torsion terms. Default: `True`.

#### `inputkeyfile`
External `.key` file whose contents are appended to the database-matched key. Default: `None`. The file is copied to `inputkey.key` for safety.

#### `prmmodfile`
Comma-separated list of parameter-modification patch files (`.mod`). See `ParameterFiles/dma4_hfe2023.mod`. Default: empty.

#### `quickdatabasesearch`
Use a faster, less exhaustive database-matching strategy. Default: `False`.

#### `mmbondtol`
Percent tolerance on MM-vs-QM bond-length deviation (`|qm-mm|*100/qm`) before flagging. Default: `0.5` (%).

#### `mmangletol`
Percent tolerance on MM-vs-QM equilibrium-angle deviation (`|qm-mm|*100/qm`) before flagging. Default: `0.5` (%).


### 10. Torsion Parameterization

#### `dontdotor`
Skip torsion parameterization entirely. Default: `False`.

#### `dontdotorfit`
Run torsion QM but skip the actual fitting. Default: `False`.

#### `rotalltors`
Force a QM scan of every rotatable bond, even those matched in the database. Default: `False`.

#### `toroptmethod`
Method for torsion-scan optimizations. Values: `xtb`, `ANI`, `AMOEBA`, or any QM method recognized by Psi4/Gaussian. Default: `xtb`.

#### `toroptmethodlist`
Comma-separated list of opt methods to compare torsion profiles across. Default: `[toroptmethod]`.

#### `torspmethod`
Method for torsion single-point energies. Default: `wB97X-D`.

#### `torspmethodlist`
Comma-separated list of SP methods to compare. Default: `[torspmethod]`.

#### `toroptbasisset`
Basis set for torsion optimizations. Default: `6-311G`.

#### `torspbasisset`
Basis set for torsion single points. Default: `6-311+G*`.

#### `toroptpcm`
PCM for torsion opt: `True` / `False` / `auto`. Default: `auto` (`-1`).

#### `torsppcm`
PCM for torsion SP: `True` / `False` / `auto`. Default: `auto` (`-1`).

#### `xtbmethod`
GFN-xTB Hamiltonian version (1, 2, or FF). Default: `2` (GFN2-xTB).

#### `xtbtorresconstant`
xTB torsion-restraint force constant (kcal/mol/deg²). Default: `5`.

#### `torsionrestraint`
Force constant (kcal/mol/rad²) restraining non-scanned dihedrals around the active rotatable bond. Default: `1641.4` (0.5 kcal/mol/deg²).

#### `torsionprmrestraintfactorL1`
L1 regularization weight on torsion-parameter magnitudes during fitting. Default: `0.1`.

#### `torsionprmrestraintfactorL2`
L2 regularization weight on torsion-parameter magnitudes during fitting. Default: `0.1`.

#### `foldnum`
Number of Fourier terms used in the torsion fit. Default: `3`.

#### `tordatapointsnum`
Number of grid points sampled along each dihedral. Default: `12` (auto-increases if more torsion parameters need fitting).

#### `defaultmaxtorsiongridpoints`
Upper bound on torsion-scan grid points. Default: `40`.

#### `onlyrotbndslist`
Restrict torsion sampling to listed bonds. Format: `i j, k l, ...`. Default: empty.

#### `dontrotbndslist`
Skip torsion sampling on listed bonds (useful for AA backbone). Format: `i j, k l, ...`. Default: empty.

#### `onlyrottortorlist`
Restrict tortor sampling to listed triples (`b c d` from `a b c d e`). Default: empty.

#### `onlyfittorstogether`
Comma-separated tuples of torsion indices to fit jointly. Default: empty.

#### `refinenonaroringtors`
Pucker non-aromatic rings and refine the in-ring torsion parameters. Default: `False`.

#### `nonaroringtor1Dscan`
Cut and 1D-scan a non-aromatic-ring torsion (instead of transferring alkane parameters). Default: `False`.

#### `fitfirsttorsionfoldphase`
Allow a phase offset on the first cosine term of the torsion fit. Default: `False`.

#### `maxtorresnitrogen`
Number of dihedral restraints per N-containing rotatable bond. Default: `2`.

#### `tortor`
Sample the 2D torsion-torsion potential-energy surface for adjacent rotatable dihedrals. Matched only when `only` is not in the line so it does not collide with `onlyrottortorlist`. Default: `False`.

#### `boltzmantemp`
Effective temperature (in kcal/mol) used for Boltzmann weighting `exp(-E/T)` of torsion conformers during fitting. Default: `8`.

#### `suppresstorfiterr`
Demote torsion-fit errors to warnings. Default: `False`.

#### `maxRMSD`
Maximum allowed RMSD (Å) between QM and MM optimized geometries. Default: `1`.

#### `maxRMSPD`
Maximum allowed RMS potential deviation in the ESP fit. Default: `1`.

#### `maxtorRMSPD`
Maximum allowed RMS deviation between QM and MM torsion energy profiles. Default: `1`.

#### `maxtorRMSPDRel`
Relative version of `maxtorRMSPD`. Default: `0.1`.


### 11. Fragmenter

#### `dontfrag`
Disable fragmenting; parameterize the whole molecule as-is. Default: `False`.

#### `maxgrowthcycles`
Maximum fragment-growth iterations before stopping. Default: `4`.

#### `WBOtol`
Maximum allowed Wiberg-bond-order difference between parent and fragment. Default: `0.05`.

#### `smallmoleculefragmenter`
Use the small-molecule-tuned fragmentation path. Default: `False`.

#### `skipespfiterror`
Continue past ESP-fit errors raised in fragment jobs. Default: `False`.

#### `fennixmodelname`
Name of the Fennix neural-network model checkpoint used for fragment energies (resolved as `<fennixmodeldir>/<fennixmodelname>.fnx`). Default: `ani2x_model0`.


### 12. Van der Waals Parameterization

#### `dovdwscan`
Refine database vdW parameters by fitting to ab initio dimer scan data. Default: `False`.

#### `onlyvdwatomlist`
Restrict vdW refinement to listed atom indices. Default: `None`.

#### `onlyvdwatomindex`
Refine a single vdW atom (by 1-indexed atom number). Default: `None`.

#### `fixvdwtyperadii`
Comma-separated vdW types whose radii are frozen during fitting. Default: empty.

#### `vdwprmtypestofit`
Which vdW parameters to fit: subset of `S` (sigma/radius), `T` (well depth), `R` (reduction factor). Default: `['S', 'T']`.

#### `fitred`
Also fit the hydrogen reduction factor. Default: `False`.

#### `addlonepairvdwsites`
Add lone-pair vdW sites for halogens/heteroatoms with sigma holes. Default: `False`.

#### `accuratevdwsp`
Use MP2/CBS counterpoise on aug-cc-pV[TQ]Z for dimer SPs (slow but accurate). Default: `False`.

#### `homodimers`
Generate homodimer scans (same molecule + same molecule). Default: `False`.


### 13. QM Output File Names

These override default file names for resuming or feeding pre-computed QM output. All default to internally generated paths derived from the molecule name (see `assign_filenames` in `poltype.py`).

#### `optlog`
Geometry-optimization log file.

#### `dmalog`
DMA single-point log file.

#### `esplog`
ESP single-point log file.

#### `dmafck`
DMA Gaussian formatted-checkpoint (fchk) file.

#### `espfck`
ESP fchk file.

#### `dmachk`
DMA chk file.

#### `espchk`
ESP chk file.

#### `formchk`
Formatted-checkpoint output file produced by `formchk`.

#### `gdmaout`
GDMA output file containing distributed multipoles.


### 14. Tolerances and Error Reporting

Most tolerance keywords are documented in their parent sections above. They are listed here together for cross-reference:

- `chargethreshold`, `dipolethreshold`, `quadrupolethreshold` — magnitude thresholds for flagging "big multipole" atoms.
- `dipoletol`, `absdipoletol` — molecular-dipole tolerances (both relative and absolute gates must trigger).
- `mmbondtol`, `mmangletol` — MM-vs-QM bonded-geometry percent tolerances.
- `maxRMSD`, `maxRMSPD`, `maxtorRMSPD`, `maxtorRMSPDRel` — geometry / potential / torsion RMS limits.
- `WBOtol` — Wiberg bond-order tolerance for fragments.
- `suppressdipoleerr`, `suppresstorfiterr`, `skipespfiterror` — demote specific error classes to warnings.


### 15. Legacy / External-Workflow Keywords

These keywords are parsed by `poltype.py` but are not used by the core parameterization pipeline. They are reserved for downstream solvation / BAR free-energy workflows and ligand-binding pipelines that consume Poltype output.

#### `uncomplexedproteinpdbname`
PDB filename of the apo protein structure for binding-affinity setup.

#### `listofligands`
Comma-separated list of ligand identifiers belonging to a binding study.

#### `templateligandfilename`
Template ligand structure (SDF/MOL) for alignment in a binding setup.

#### `templateligandxyzfilename`
Template ligand Tinker XYZ used for atom-index mapping.

#### `barinterval`
Lambda-window interval for Bennett Acceptance Ratio (BAR) free-energy calculations.

#### `targetenthalpyerror`
Convergence threshold for enthalpy error in liquid-property fitting.

#### `targetdensityerror`
Convergence threshold for density error in liquid-property fitting.

#### `enthalpyrelativeweight`
Weight of the enthalpy term in liquid-property fitting.

#### `densityrelativeweight`
Weight of the density term in liquid-property fitting.
