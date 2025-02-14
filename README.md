# Automation of AMOEBA Polarizable Force Field for Small Molecules - Poltype 2

## Obejective
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

[Parameterization Input Preparation](#parameterization-input-preparation)

[Ligand Protonation State Generation](#ligand-protonation-state-generation)

[Minimum Example Usage Parameterization](#minimum-example-usage-parameterization)

[Recommended-Poltype-inputs](#recommended-poltype-inputs)

[Default Resource Consumption](#default-resource-consumption)

[💻 Advanced Program Usage](README/README_HELP.MD)

[Parameterization Output Files](README/README_OUTPUT.MD)

* [Final XYZ File](README/README_OUTPUT.MD)

* [Final Key File](README/README_OUTPUT.MD)

* [Poltype Log File](README/README_OUTPUT.MD)

* [OPENME Plots](README/README_OUTPUT.MD)

[Parameterization Sanity Checks](#parameterization-sanity-checks)

[Parameterization Examples](Examples/Parameterization)

[Automated AMOEBA Ligand Parameterization How It Works](README/README_PROTOCOL.MD)


[💻 Advanced Program Usage](README/README_HELP.MD)

[Automated AMOEBA Molecular Dynamics and Free Energy Prediciton How It Works](README/README_PROTOCOL.MD)



---------------------------------------------------------------------------------------------


    
### Parameterization Input Preparation
* The input structure can be given as an filetype with bond order information (sdf,mol,mol2,etc..)
* Cartesian XYZ, tinker XYZ and PDB is not recommended, this will be converted to an SDF file and bond orders will be guessed based on atomic distances.
* If 2D structure is given, poltype will generate 3D coordinates for you.
* Total charge is determined by computing the formal charge for each atom.
* Formal atom charge will be assigned via the input number of bonds and bond order for surrounding bonds and element of each atom. 
* Optional keywords exist to add missing hydrogens.
* If carbon or nitrogen are negatively charged, then Poltype will assume you meant to have hydrogens on those atoms and protonate to neutral charge. For all other elements, if there exists a formal charge (due to input bond order and valence electrons of element), then poltype will detect the formal charge for you and compute total charge from sum of each formal atom charge. Warning message is printed when hydrogen atoms are added. 
* Special radical charge states require additional information in the input file specifying which atom is a radical. 
* If you do not want charged atoms to be protonated then use ``addhydrogentocharged=False``
### Ligand Protonation State Generation
* Dominant ionization states at pH 7 are enumerated and SDF files are generated via Dimorphite-DL (IonizationState_0.sdf,IonizationState_1.sdf,..). 
* Tautomer states are enumerated via rdkit and the first tautomer is the canonical tautomer TautomerState_0.sdf
* Use ``genprotstatesonly`` to quit program after generating dominant ionization states at pH=7 and tautomers.

### Minimum Example Usage Parameterization

__All input arguments are specified in poltype.ini file__
```
structure=methylamine.sdf
```
* Navigate to directory containing poltype.ini and .sdf file, and run:

```shell
nohup python /path_to_poltype/poltype.py &
```

```final.xyz``` and ```final.key``` are the resulting structure and parameter files you will need.
* After poltype finishes, check the ``OPENME`` folder for torsion fitting and ESP fitting results. 

### Recommended Poltype inputs

By default, Poltype uses the old GDMA algorithm with non-bonded parameters mostly from amoeba09.

To use the new GDMA algorithm and non-bonded parameters adjusted to reproduce hydration free energies:
```
new_gdma=True
gdmacommand_Radius_S=0.80
prmmodfile=dma4_hfe2023
```
[Adjustments inlcude multipole and vdw scaling](ParameterFiles/dma4_hfe2023.mod)

To use the new GDMA algorithm with default GDMA radius, no scaling of multipole or vdW parameters:
```
new_gdma=True
```

### Default Resource Consumption
* By default, Poltype computes the number of fragment poltype jobs (or any QM job if fragmenter is not being used) to run in parallel as the floor function of the input number of cores divided by the number of cores per job (default of 2). 
* RAM, cores, and disk space can all be detected and a consumption ratio of 80% is used by default. Fragment RAM, disk, and cores are divided evenly by the number of Poltype jobs in parallel. 




### Parameterization Sanity Checks
* MM = Molecular Mechanics (AMOEBA model), QM = Quantum Mechanics
* Check for 2D coordinates and generates 3D coordinates at begining of program 
* Check to ensure refined multipoles enable MM to model the QM potential grid well and raises error if not
* Check to ensure QM and MM dipoles are very similar and raise error if not
* Check to ensure the final minimized MM structure is similar to the geometry optimized QM structure and raises error if not
* Check for any missing van der Waals at end of program parameters and raises error
* Check for any missing multipole parameters at end of program and raises error
* Check for any zeroed out torsion parameters in final key file. Checks if fragmenter is transferring torsion properly.

