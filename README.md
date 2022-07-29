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
* Host-guest Modelling Features
    * Docking with GOLD, AutoDock4, AutoDock Vina, Vinardo
    * Protein-ligand interaction profiler and visualization with BINANA and ProLif
    * Missing residues & loop modeling with Modeller
    * Protonation assignment via propka/pdb2pqr
    * Tinker box set up
    * BAR Free energy estimation file setup



[💻 Program Installation](README/README_INSTALL.MD)


### Automated AMOEBA Ligand Parameterization 

[Parameterization Input Preparation](#parameterization-input-preparation)

[Ligand Protonation State Generation](#ligand-protonation-state-generation)

[Minimum Example Usage Parameterization](#minimum-example-usage-parameterization)

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

### Automated AMOEBA Molecular Dynamics and Free Energy Prediciton

<img src="https://florentbarbault.files.wordpress.com/2010/09/mova1.gif" width="20%">

* [1](#citations)

[Protein Ligand Interaction Visualization](#protein-ligand-interaction-visualization)

[Protein Input Preparation](#protein-input-preparation)

[Minimum Input Example Docking](#minimum-input-example-docking)

[Molecular Dynamics Input Preparation](#molecular-dynamics-input-preparation)

[Minimum Input Example Binding Free Energy](#minimum-input-example-binding-free-energy)

[Minimum Input Example Solvation Free Energy](#minimum-input-example-solvation-free-energy)

[Minimum Input Example Neat Liquid Simulation](#minimum-input-example-neat-liquid-simulation)

[Molecular Dynamics Sanity Checks](#molecular-dynamics-sanity-checks)

[Free Energy Output Files](#free-energy-output-files)

*  [Binding Free Energy Table Output Example](#binding-free-energy-table-output-example)

*  [Hydration Free Energy Table Output Example](#hydration-free-energy-table-output-example)

[Hydration Free Energy Examples](Examples/HydrationFreeEnergy)
 
[Binding Free Energy Examples](Examples/BindingFreeEnergy)

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



### Protein Ligand Interaction Visualization

<img src="README/Images/3DProLig.png" width="50%">
<img src="README/Images/2DProLig.png" width="50%">


* Navigate to folder VisualizationNotebooks
* Be sure to install the conda environment (notebookenvironment.yml),``conda env create -f notebookenvironment.yml``. 
* Activate the conda environment ``conda activate pymolenv``
* Move protein-ligand complexed PDB to VisualizationNotebooks folder.
* Ensure that the residue labels for ligand are labeled as ``LIG``.
* Launch the jupyter notebook ``jupyter-notebook Protein-Ligand-Interactions.ipynb``
* Input your complexed PDB name into variable ``ligandreceptorfilename``
* Run the cell and the output from BINANA interaction profiler and ProLif is shown. 


### Protein Input Preparation
* If your protein PDB has missing residues, poltype wraps Modeller to fill in missing residues for you (experimental). Need to use keyword ``pdbcode``. This will download the PDB for you, and call modeller, then quit the program after PDB has been filled in and optimized with Modeller. 
* An academic license file is required after installation of modeller ``conda install -c salilab modeller``, the screen will prompt you which file to insert your license key in.
* After this check results of output PDB.
```
pdbcode=5l2t
```
* Remove the keyword ``pdbcode`` from poltype input file, if you wish to perform further computations. 
* Protonation state assingment and adding ligand to the protein pocket are next steps for computing binding simulations.
* ``usepdb2pqr`` keyword can be used with ``uncomplexedproteinpdbname`` to estimate pKa values of titratable residues via propka and then protonate the PDB for you. The output PDB will have extension "_final.pdb".
* pdb2pqr can be installed via ``conda install -c conda-forge pdb2pqr`` or use yaml file
```
uncomplexedproteinpdbname=5l2t_filled.BL00020001.pdb
```
* After adding any missing residues and assigning the protonation state, the ligand needs to be added to the protein pocket and given as input ``complexedproteinpdbname`` for binding computations.


### Minimum Input Example Docking
* If you want to use GOLD, ensure that the bin folder is in your PATH
* Complexed protein with ligand needs to be provided (this should already be protonated before providing as input)
* For generating PDBQT files with AutoDock4 and AutoDock Vina, a seperate python 2 environment (dockingprep.yml) needs to be installed.
* Optional keywords exist to change docking grid center (default is center of ligand from the input protein-ligand complex), docking grid length (how far grid extends), grid spacing (controls point density on grid) and the number of poses generated (default 10).  
* Final scores, structures and rankings are given in a file called DockingReport.txt 

```
complexedproteinpdbname=complex.pdb
listofligands=lig1.sdf,lig2.sdf
usead4=True
usevina=True
usevinardo=True
usegold=True
```



### Molecular Dynamics Input Preparation
* If Tinker9 executables are in PATH, then program will switch to using analyze9,dynamic9,minimize9, the GPU executables.
* Make a seperate folder from where parameterization files from poltype were made (with new poltype.ini file too)
* Ligand XYZ and key files are required (such as final.xyz and final.key from Poltype parameterization).
* Besides the ``ligandxyzfilenamelist`` for all ligands, there is also ``annihilateligandxyzfilenamelist`` which tells the program which ligand types to annihilate. By default, annihilateligandxyzfilenamelist=ligandxyzfilenamelist.
* If there are duplicate molecules (same parameters and XYZ), make copies of each XYZ and key and provide as inputs. The program will then change the type numbers so as none of the individual molecules have overlapping type numbers. This way can choose to disappear one or more of the ligands with the same parameters (but different types). If there are duplicate ligands in your system, please make sure the list order of ``ligandxyzfilenamelist`` and ``annihilateligandxyzfilenamelist`` are the same order as ligands that occur in the input PDB file. This way the program can determine which of the duplicates in ``annihilateligandxyzfilenamelist`` to disappear and distinguish those from other duplicates in ``ligandxyzfilenamelist`` that are not in ``annihilateligandxyzfilenamelist`` 
* Waters and ions within 8 angstroms of any atom in ligand - protein PDB will automatically be detected and added to protein-ligand XYZ file.
* For binding free energy compuations, either a host PDB or premade tinker XYZ is required 
* Inputs are inside poltype.ini
* Make sure pdb files (complexed and uncomplexed) have no missing residues or atoms.
* Charge is read from input XYZ files generated.
* Make sure if using custom receptor parameters, then either adding to keyfilenamelist or in prmfilepath
* Use submitlocally=False if you do not wish to submit dynamics jobs locally. By default this is already False for production dynamics and for BAR, for minimzation and equilbriation, by default this is True and jobs are submited locally. Then program will wait for you to complete the jobs in text file (such as _proddynamicsjobs.txt).
* For HFE, if your ligand is charged and you want to compute the salt hydration free energy, add "salthfe=True"
* Some files containing .out extensions are used to determine if dynamics has already been done to prevent rerunning dynamics each time program is called. Tinker keys files are always recreated everytime program is rerun. 
* Minimized and equilibriated PDB (ligand complexed with host) and XYZ structures for viewing in pymol will have extension _pymol in name. Also PDB trajectory files are generated for complexed host-ligand (if input complexedproteinpdbname is given) for equilbriated ARC and production dynamics ARC with full electrostatics and vdW interactions. PDB trajectory files are in same location as ARC file. 
* If your dynamic jobs are killed prematurely, upon rerunning poltype it will read the number of steps in your output dynamics files and determine the number of new steps needed to be taken to finish dynamics.

#### Minimum Input Example Binding Free Energy

```
complexedproteinpdbname=complex.pdb 
keyfilenamelist=Mol1.key , Mol2.key
ligandxyzfilenamelist=Mol1.xyz , Mol2.xyz
```
or

```
receptorligandxyzfilename=complex.xyz
prmfilepath=prmfile # with absolute path for receptor
keyfilenamelist=complex.key 
ligandxyzfilenamelist=ligand.xyz
```

#### Minimum Input Example Solvation Free Energy

```
keyfilenamelist=Mol.key
ligandxyzfilenamelist=Mol.xyz
```

#### Minimum Input Example Neat Liquid Simulation
```
density=997
keyfilenamelist=Mol.key
ligandxyzfilenamelist=Mol.xyz
equilibriatescheme=50,100,150,200,300,300
```

* Navigate to directory containing poltype.ini, and run:

```shell
nohup python /path_to_poltype/poltype.py &
```


### Molecular Dynamics Sanity Checks
* Check for missing parameters
* Check total charge of all boxes for binding free energy alchemical perturbation have a net zero charge

### Free Energy Output Files


#### Binding Free Energy Table Output Example
* ΔGˢᵒˡᵛ = Change in solvation free energy (does not contain gas phase component by default)
* ΔGˢᵒˡᵛᵉʳʳ = Change in solvation free energy error
* ΔGᶜᵒᵐᵖᶜᵒʳʳ = Change in complexation free energy with analytical correction
* ΔGᶜᵒᵐᵖᵘⁿᶜᵒʳʳ = Change in complexation free energy with out analytical correction
* ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ = Change in complexation free energy error with analytical correction 
* ΔGᵃⁿᵃᶜᵒᵐᵖᶜᵒʳʳ = Analytical correction to complexation free energy
* ΔGᵇᶦⁿᵈᶜᵒʳʳ = Change in binding free energy corrected
* ΔGᵇᶦⁿᵈᶜᵒʳʳᵉʳʳ = Change in binding free energy corrected error
* ΔGˢᵒˡᵛᵉˡᵉ = Change in solvation free energy electrostatic component
* ΔGˢᵒˡᵛᵛᵈʷ =  Change in solvation free energy van der Waals component
* ΔGᶜᵒᵐᵖᵉˡᵉ = Change in complexation free energy electrostatic component
* ΔGᶜᵒᵐᵖᵛᵈʷ =  Change in complexation free energy van der Waals component
* ΔHˢᵒˡᵛ = Change in solvation enthalpy 
* ΔHˢᵒˡᵛᵉʳʳ = Change in solvation entropy error
* ΔSˢᵒˡᵛ = Change in solvation entropy
* ΔSˢᵒˡᵛᵉʳʳ= Change in solvation entropy error
* ΔHᶜᵒᵐᵖ = Change in complexation enthalpy
* ΔHᶜᵒᵐᵖᵉʳʳ = Change in complexation enthalpy error
* ΔSᶜᵒᵐᵖ = Change in complexation entropy
* ΔSᶜᵒᵐᵖᵉʳʳ = Change in complexation entropy error
* ΔHᵇᶦⁿᵈ = Change in binding enthalpy
* ΔHᵇᶦⁿᵈᵉʳʳ = Change in binding enthalpy error
* ΔSᵇᶦⁿᵈ = Change in binding entropy
* ΔSᵇᶦⁿᵈᵉʳʳ = Change in binding entropy error

<img src="README/Images/GibbsBindingTable.png" width="100%">

* Gibbs free energy table 
* Gibbs_Free_Energy_Change_Table.csv

<img src="README/Images/Enthalpy_Entropy.png" width="100%">	

* Enthalpy, Entropy and Gibbs free energy table
* Enthalpy,_Entropy,_Gibbs_Energy_Change_Table.csv

<img src="README/Images/SolvationSimTable.png" width="100%">

* Generic simulation information
* Solvation_Simulation_Info_Table.csv

<img src="README/Images/ComplexationSimTable.png" width="100%">

* Generic simulation information
* Complexation_Simulation_Info_Table.csv



#### Hydration Free Energy Table Output Example
* ΔGˢᵒˡᵛ = Change in solvation free energy
* ΔGˢᵒˡᵛᵉʳʳ = Change in solvation free energy error
* ΔHˢᵒˡᵛ = Change in solvation enthalpy 
* ΔHˢᵒˡᵛᵉʳʳ = Change in solvation enthalpy error
* ΔSˢᵒˡᵛ = Change in solvation entropy
* ΔSˢᵒˡᵛᵉʳʳ = Change in solvation entropy error
* ΔGˢᵒˡᵛᵉˡᵉ = Change in solvation free energy electrostatic component
* ΔGˢᵒˡᵛᵛᵈʷ =  Change in solvation free energy van der Waals component
* ΔGᵉˡᵉᵍᵃˢ =  Change in solvation free energy electrostatic component gas phase
* ΔGᵉˡᵉˢᵒˡ =  Change in solvation free energy electrostatic component solution phase
* ΔGᵛᵈʷᵍᵃˢ = Change in solvation free energy van der Waals component gas phase
* ΔGᵛᵈʷˢᵒˡ = Change in solvation free energy van der Waals component solution phase
* ΔGˢᵒˡ = Change in solvation free energy solution phase
* ΔGᵍᵃˢ = Change in solvation free energy gas phase
* ΔGˢᵒˡᵛᶠʷᵈ = Change in forward solvation free energy (A->B)
* ΔGˢᵒˡᵛᵇʷᵈ =  Change in backword solvation free energy (B->A)
* SolvOverlap = Solvation overlap between neigboring states, value 0-1





<img src="README/Images/GibbsHFE.PNG" width="100%">

* Gibbs free energy table 
* Gibbs_Free_Energy_Change_Table.csv

<img src="README/Images/EnergyHFE.PNG" width="100%">	

* Enthalpy, Entropy and Gibbs free energy table
* Enthalpy,_Entropy,_Gibbs_Energy_Change_Table.csv

<img src="README/Images/BoxInfoHFE.PNG" width="100%">

* Generic simulation information
* Solvation_Simulation_Info_Table.csv

<img src="README/Images/BarResultsHFE.PNG" width="100%">	

* Individual BAR step free energy computations
* Useful for troubleshooting if free energy is incorrect
* BARResults.csv






