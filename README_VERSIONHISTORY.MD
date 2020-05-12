## Version History

[2.0.0](https://github.com/pren/poltype/tree/346ec870775c15e7279238ce8c9bba8bd267e36e)
* Poltype is modularized for code organization.
* Unittesting is available for maintaining code robustness.
* Updated documentation, Examples folder is included.
* Various high symmetry molecules are handles, see Example folder.
* Torsion spinning is parallelizable.
* Modified Amino Acids can be swapped with canonical amino acids in protein PDBs and parameters can be generated for side chain then automatically stitched together with protein forcefield.

[2.0.1](https://github.com/pren/poltype/tree/3bc6038be63c2e91a06f06b24f96a8e0fa462b91)
* Fix bug in scratch directory assignment for Gaussian

[2.0.2](https://github.com/pren/poltype/tree/2ad7ee5a7cb4c8185e96e6299c3b3db9d97fcb4f)
* rdkit does not handle sanitization of some aromatic compounds during torsion spining. This has been turned off.

[2.0.3](https://github.com/pren/poltype/tree/1e24c84bcc98938af52a5fe834a73d13281def2e)
* Fix bug that prevented post fitting of multiple torsion in molecule (unittesting only tested for small molecules ~ 1 torsion).

[2.0.4](https://github.com/pren/poltype/tree/eb3bfd4209d70d8931d4cc88fddb02c2098921b0)
* Now when quantum log file exists but is blank, POLTYPE will skip this data point from torsion fitting
* Previously it was assumed if optimization log file exists, then the rest of quantum (SP) finished as well. This assumption has been removed
* function for generating SP for torsion was referecing espmethod and espbasisset instead of torspmethod and torspbasisset

[2.0.5](https://github.com/pren/poltype/tree/4d640adbf25a7cd1e060546851f592937681a080)
* If RMSPD.txt already existed then this file was not recreated, this is an issue if the first time created, tinker crashed. This assumption has been removed.

[2.0.6](https://github.com/pren/poltype/tree/51684891f5993d17b16302b58a52022bc37510cc)
* Now when Tinker version number is not correct a ValueError will be raised
* Now if an external command (for example calling Gaussian, Psi4 etc..) fails, poltype will crash with an error instead of attempting to continue
* Changing default optmethod and toropmethod to wB97X-D
* fragmenter not working yet, add dontfrag to poltype.ini files, need to push updates though to keep poltype.py consistent with other changes

[2.0.7](https://github.com/pren/poltype/tree/734362e8baa560761a650611381f427a30789753)
* Default SP for torsion is set back to MP2 for accuracy

[2.0.8](https://github.com/pren/poltype/tree/ca98ccf914cc8570b833460f9ad79ed051ec1c80)
* Fragmenter turned off by default
* anglep not assigned unless middle atom is connected to three atoms to be consistent with TINKER
* Using relative error for QM and MM Dipole comparison since the MMDipole scales worse for larger molecules
* Fix bug related to making QM logs blank files and skipping them

[2.0.9](https://github.com/pren/poltype/tree/218752f6cf63771f0093140627d7da1bc48ef43a)
* If bond is a ring bond or at least one atom in a bond is a ring atom, cannot spin dihedral angle according to rdkit, so just transfer from database
* Bug fix with charged molecules, solvate GK was turned on during ESP fitting, this is turned on only during torsion fitting

[2.1.0](https://github.com/pren/poltype/tree/2cbfb95ed2b1224e6b980ebbce00a3d5e870f756)
* If running torsion serially on local host instead of crashing poltype when Gaussian/Psi4 job fails, just continue and attempt to do torsion fitting without that data point

[2.1.1](https://github.com/pren/poltype/tree/2ca87c59780b542e8f8031f1d9c6a38d22180287)
* Added TARGET-DIPOLE dipole from ESP QM data to better reproduce MM QM dipole agreement. Gaussian uses center of nuclear charge frame which we convert to center of mass coordinates (Tinker is in COM Frame). Psi4 default is COM frame.

[2.1.2](https://github.com/pren/poltype/tree/93afbeb6d59951f3e47f35da3d20daa35fb85729)
* Bug fix for TARGET-DIPOLE keyword being used before ESP QM generated

[2.1.3](https://github.com/pren/poltype/tree/90ea499858f1c505bc2f4905e78a7c08a46f174e)
* Added SOLUTE keyword for GK and DD-COSMO electrostatic radii

[2.1.4](https://github.com/pren/poltype/tree/89d30cfd9d10d01713a57dfec4038b63e79a8478)
* Check RMSPD root mean square potential difference and raise error before checking dipoles, (previously supposed to do this but there was a bug)

[2.1.5](https://github.com/pren/poltype/tree/4ead5cdf6a78d45d67077976703cb94831a7d7e4)
* When extracting dipole from gaussian dipole is flipped in wrong direction, this is fixed

[2.1.6](https://github.com/pren/poltype/tree/214755fe8021d647a8a5c478d90dbe903b26879d)
* Error with generating opt.com file for modifiedresidue library with first poltype run, this has been fixed. Now when poltype finishes, ttt.xyz_2 is deleted so if rerunning ttt.xyz_3 ... isn't created, also .chk files are deleted if job is successful.

[2.1.7](https://github.com/pren/poltype/tree/6b8b89e52b0b3153b42bce3c393934d3c84f10dc)
* Torsion parameters are transferred from database if not in rotatable bond list (bug fix from previous edits of rotatable bond definition and transferable torsion in database)

[2.1.8](https://github.com/pren/poltype/tree/4fc6544b81939e4dda8662cebbe97cb2a7674837)
* Bug fix, solvate GK term removed from ttt.key for charged molecules to compute dipoles (ESP fitting done in gas phase, quantum done in gas phase)

[2.1.9](https://github.com/pren/poltype/tree/564cbb31ff3fe44bcb091b3ca71cd30a90f50898)
* Input structure for mutated residues can have hydrogens on back bone now. Bug for finding all transferable torsions is fixed. See Example ModifiedAminoAcid
* Instead of using babel to predict bonds from PDB via distances a library with SMART strings for amino acids and nucleic acids (residue_connect.txt) is used to match to PDB guess via distances but only keep matched bonds and also bonds connecting residues. The babel guess will sometimes but bonds between sidechains if they are too close to each other in the PDB file so a library was needed instead.

[2.2.0](https://github.com/pren/poltype/tree/42878c15b60011cbc3c40f7ca45c7c078aaa19ef)
* PSI4 error for torsion SP quantum jobs computing Wiberg Bond Index, forgot to pass wavefunction to oeprop class

[2.2.1](https://github.com/pren/poltype/tree/58e5b4e64fc68497e09493f5e77c57fd090d7739)
* Fix error in -pedit.txt where a space was delayed before being written out to text file causing space to be in wrong place, fixed by flushing system buffer

[2.2.2](https://github.com/pren/poltype/tree/2d44a898b0b0774e26fcc84c29ccca523944b7ce)
* Updates to README files as well as fragmenter, torsion libraries, unittesting library.

[2.2.3](https://github.com/pren/poltype/tree/bb4ea928bb9f759da42d040e50447649fc1ee4bd)
* Bug in TINKER 8.7 POLEDIT.X, POLTYPE now supports TINKER 8.2 and TINKER >=8.7

[2.2.4](https://github.com/pren/poltype/tree/8e5ba282ebe4b52716e8afc58770eb83a2f53070)
* Bug is okay now, continue using TINKER >=8.7

[2.2.5](https://github.com/pren/poltype/tree/d07670aca3543ce9fd86e0c2e4bfd37e439be521)
* multipole scaling turned off for alcohols

[2.2.6](https://github.com/pren/poltype/tree/78132c7b51a58b711d5abe290f7a980065d4dde7)
* fixed symmetry frames bug, added new Symmetry molecules cases to Examples folder, see README.MANIFEST

[2.2.7](https://github.com/pren/poltype/tree/c7092f9667bc2392019308da26e1b7bca941f692)
* Torsion QM jobs parallelized for all torsion as opposed to for a given torsion previously. Example daemon is in Examples folder, see README.MANIFEST.

[2.2.8](https://github.com/pren/poltype/tree/c1a4ba264072ed2b1a558b1543b7e889575c4539)
* Fix bug in parallel torsion

[2.2.9](https://github.com/pren/poltype/tree/75733a749041b98bffca521f78d24bf201787b2a)
* First working version of torsion fragmenter, this is turned on by default. All optimization have maxcycles of 5, since bond and angle degrees of freedom converge more rapidly than conformation. MMFF94 is used in openbabel to minimize geometry before QM optimization is done.

[2.3.0](https://github.com/pren/poltype/tree/1370038622322622fedd4b68f13e6e425c59f5e9)
* Bug with no torsion to fit but entering fit function causing crash is fixed
