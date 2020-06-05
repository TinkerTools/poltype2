* [ ] Fix poltype GUI
* [ ] Add graph visualization of function calls in poltype
* [ ] Add RaiseValue error for when child poltype job has not completed (doesnt have key_5), otherwise people use fragmenter and wont read the parent poltype log that contains all of the errors from any child jobs
* [ ] For some reason old poltype did not restrain hydrogen torsion on QM (even though MM does) , make sure it is apples to apples
* [ ] Fix unittester.
* [ ] Instal theme on GitHub page
* [ ] Add README_TROUBLESHOOTING page (copy content from README_HELP). 
* [ ] Before converting fragment to SMARTS/SMIRKS for database, remove any added hydrogens from cut bonds
* [ ] Add documentation about polarization groups.
* [ ] Docker contained for Linux OS
* [ ] Allow option to allow users to have correct formal charges assigned in input file and turn off automatic formal charge assignment in fragmenter lib.
* [ ] If freq option is turned on, need to make sure optimization doesnt stop after default 5 optmaxcycles (need opt to converge)
* [ ] Use unittester on fragmenter with different WBO tolerances to determine what tolerance is acceptable for a given tolerance in RMSD between parent QM-MM and fragment QM-MM. 
* [ ] People seem to change optmethod a lot without deleting XYZ files or .cube, fchk files causing RMSPD error. Just auto delete XYZ files (quick to remake). 
* [ ] Give input option to add hydrogen.
* [ ] Make sure gaussian deletes scratch files when job is finished (Psi4 does this).
* [ ] Make sure of key_3 already exist, dont redo fitting of esp.
* [ ] Fix bug where WithIndex and WithoutIndex.png have no index. Make sure all molecules in image either have indexes or dont. 
* [ ] Highlight the atoms that are different from molecule in previoius picture in images with many molecule.
* [ ] Add math to markdown READMEs (WBO equation)
* [ ] Add arrows to FragmentGrowth images to make it more like flow diagram.
* [ ] Remove any file operations done in bash and replace with python language (users different OS)
* [ ] Add 10% of QM-MM for torsion max_amp restrained fitting
* [ ] Remove redundant torsion parameter SMART string parameters in current database (valence).
* [ ] Raise ValueError if final key file has any torsion zeroed out (possible that database had no matches)
* [ ] For any hydrogen torsion assign prior SMART string and parameters to reduce the amount of parameters needed to fit for rotatable bonds
* [ ] Add map from Tinker Type Number to SMART string chemical enviornment
* [ ] Using the above map, autogenerate torsion, vdw ... parameters from .prm file inside poltype. Keep user added parameters in seperate dictionaries. 
* [ ] For a given node, RAM, processor number and disk usage, can we parralelize QM jobs for a single node (as opposed to just submitting one job per node?)
* [ ] For adding vdw parameters need SMART string to be certain "size" to be transferable in database (just do for user in auto vdw fitter).
* [ ] For any autofitting and fragment SMARTS generation,also include metrics for RMSED for fitting in comments (can be automated) so future users know how bad fit is 
* [ ] Fix bugs in modifedresidue.py, also fix bug giving excess charge (small) when starting with neutral input molecule.
* [ ] Allow charged fragment molecules for fragmenter when formal charge on parent detected
* [ ] Add 6 fold torsion
* [ ] add weighted torsion for fitting function
* [ ] add z-x-y multipole frame for chiral atoms
* [ ] Test fragmenter on drugbank database (dont paramterize until the rest of POLTYPE torsion is complete)
* [ ] Add multiple conformers for multipole fitting
* [ ] Allow ring torsions to be puckered and obtain QM data for torsion paramterization
* [ ] Allow 2D dihedral scans to fit a tor-tor coupling term too for torsions that are adjacent to each other. Numerical spline should be used for functional form of tor-tor.
* [ ] Allow user to specify if they want cosine terms for torsion or if they want numerical QM-MM spline (partitioned and then we put error on heaviest torsion)
* [ ] Add support for AMOEBA+ parameters
* [ ] Add automated vdw QM surface scan from water probes. Add this BEFORE torsion fitting, so that you always do torsion last. If user decides to do vdw fitting after key_5 already generated need to remove the key5 then and refit torsion.
* [ ] Add comments to every function in poltype with the fields below
    Intent: 
    Input:
    Output: 
    Referenced By: 
    Description 
* [ ] After torsion fitting refine initial guess for bond angles, etc with frequencies from valence
