* [ ] Add math to markdown READMEs
* [ ] Add TORSION_README.MD to explain torsion fitting plots and how to interpret them.
* [ ] Remove any file operations done in bash and replace with python language (users different OS)
* [ ] Allow charged fragment molecules for fragmenter when formal charge on parent detected
* [ ] Add 6 fold torsion
* [ ] add weighted torsion for fitting function
* [ ] Test fragmenter on drugbank database (dont paramterize until the rest of POLTYPE torsion is complete)
* [ ] Add multiple conformers for multipole fitting
* [ ] Allow ring torsions to be puckered and obtain QM data for torsion paramterization
* [ ] Allow 2D dihedral scans to fit a tor-tor coupling term too for torsions that are adjacent to each other. Numerical spline should be used for functional form of tor-tor.
* [ ] Add support for AMOEBA+ parameters
* [ ] Add automated vdw QM surface scan from water probes.
* [ ] Remove the need to install gdma and use gdma built in to Psi4
* [ ] Remove the need to install Tinker CPU independelty and install openmm with conda, then call tinker via openmm instead of command line options
* [ ] Once all depdencies are conda based, add poltype to anaconda to make installation very simple
