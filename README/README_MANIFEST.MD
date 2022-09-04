## Manifest

* [Poltype Files](#poltype-files)
* [Poltype Modules](#poltype-modules)
* [Poltype Examples](#poltype-examples)

### POLTYPE Files

* .sdf               Input structure file.
* -opt_1.com         Gaussian input file for first geometry optimization.
* -opt_1.psi4        Psi4 input file for first geometry optimization.
* -opt_1.log         Gaussian or Psi4 output file for first geometry optimization.
* -opt_1.chk         Unreadable electron density from first geometry optimization.
* -opt_1.fchk        Human readable electron density from first geometry optimization.
* -opt_2.com         Gaussian input file for second geometry optimization.
* -opt_2.psi4        Psi4 input file for second geometry optimization.
* -opt_2.log         Gaussian or Psi4 output file for second geometry optimization.
* -opt_2.chk         Unreadable electron density from second geometry optimization.
* -opt_2.fchk        Human readable electron density from second geometry optimization.
* -dma.com           Gaussian input file for QM SP calculation to obtain density for running distributed multipole analysis.
* -dma.psi4          Psi4 input file for QM SP calculation to obtain density for running distributed multipole analysis.
* -dma.log           Gaussian or Psi4 output file for QM SP calculation to obtain density for running distributed                         multipole analysis.
* -dma.chk           Unreadable density for running distributed multipole analysis.
* -dma.fchk          Human readable density for running distributed multipole analysis.
* -esp.com           Gaussian input file for QM SP calculation to obtain density for running electrostatic potential                      fitting.
* -esp.psi4          Psi4 input file for QM SP calculation to obtain density for running electrostatic potential fitting.
* -esp.log           Gaussian or Psi4 output file for QM SP calculation to obtain density for running electrostatic                       potetnial fitting.
* -esp.chk           Unreadable density for running electrostatic potential fitting.
* -esp.fchk          Human readable density for running electrostatic potential fitting.
* .gdmain            Input for GDMA.
* dma.punch          Intermediate file produced from GDMA.
* .gdmaout           Output from GDMA.
* .key               First keyfile containing atom indexes, multipoles from GDMA output and polarize terms.
* .grid              Textfile specifying coordinates of where electrostatic potential will be evaluated.
* .cube              Potential from QM SP density at specified grid coordinates, in format not recognized by Tinker.
* .pot               Potential from cube file at specified grid coordinates in Tinker readable format.
* -peditin.txt       Input for poledit to change frames as specified by POLTYPE.
RMSPD.txt           Textfile computing the root mean square potential difference between potential from QM and from MM.
* -groups.txt        Input for amvmpoles.pl to reduce atom typing.
* .xyz_2             Structure using symmetry types rather than atom index types.
* .key_2             Second keyfile using same parameters as *.key but with symmetry type numbers rather than atom index types.
* .key_3              Third keyfile containing multipole parameters from electrostatic potential fitting.
* .key_4              Fourth key file using same parameters as *.key_3 but also includes valence parameters from database.
* qm-torsion          Folder containing QM,MM optimized structures and SP computations to evaluate energy along dihedral angle.
* .key_5              Fifth key file using same parameters as *.key_4 but with fitted torsion parameters.
* MMDipole.txt        Contains information on MM dipole.
* QMDipole.txt        Contains information on QM dipole, only written out for Gaussian currently.
* final.xyz           Final QM optimized structure.
* final.key           Final key file.
* final.xyz_2         Final MM optimized structure from QM.
* -poltype.log        Poltype workflow, computations and external calls are shown in here.
* torsionprmguess.txt Used for storing torsion parameters needed for ring puckering of non-aromatic torsion                                refinement in rings.py module



### qm-torsion
b_c_* = the middle two atoms involved in all torsions for a given torsion scan set, (a1-b1-c1-d1,a2-b2-c2-d2,..), would translate to b_c_*=b1-c1_b2-c2

* -opt-b_c_*_angle_preQMOPTprefit.xyz           Tinker XYZ from previous Gaussian log file,before current QM optimzation                                                   and torsion fitting.

* -opt-b_c_*_angle_preQMOPTprefit_xyzformat.xyz Tinker XYZ from previous Gaussian log file,before current QM optimzation                                                   and torsion fitting.

* -opt-b_c_*_angle_preQMOPTprefit.key           Tinker key using .key_4 parameters but with additional torsion restraints                                                  to fix other torsions during minimization, before current QM optimzation                                                   and torsion fitting.

* -opt-b_c_*_angle_preQMOPTprefit.xyz_2         Minimized Tinker XYZ, before current QM optimzation and torsion fitting.

* -opt-b_c_*_angle_preQMOPTprefit.out           Output from Tinker minimization, before current QM optimzation and torsion                                                 fitting.

* -opt_1-b_c_*_angle.com                      Gaussian input for first geometry optimization.
* -opt_1-b_c_*_angle.psi4                     Psi4 input for first geometry optimization.
* -opt_1-b_c_*_angle.chk                      Unreadable density information from first geometry optimization.
* -opt_1-b_c_*_angle.log                      Gaussian or Psi4 output from first geometry optimization.

* -opt_2-b_c_*_angle.com                      Gaussian input for second geometry optimization.
* -opt_2-b_c_*_angle.psi4                     Psi4 input for second geometry optimization.
* -opt_2-b_c_*_angle.chk                      Unreadable density information from second geometry optimization.
* -opt_2-b_c_*_angle.log                      Gaussian or Psi4 output from second geometry optimization.

* -opt-b_c_*_angle_postQMOPTprefit.xyz          Tinker XYZ from previous Gaussian log file,after current QM optimzation                                                    and torsion fitting.

* -opt-b_c_*_angle_postQMOPTprefit.key          Tinker key using .key_4 parameters but with additional torsion restraints                                                  to fix other torsions during minimization, after QM optimzation and                                                        torsion fitting.

* -opt-b_c_*_angle_postQMOPTprefit.xyz_2        Minimized Tinker XYZ, after QM optimzation and torsion fitting.

* -opt-b_c_*_angle_postQMOPTprefit.out          Output from Tinker minimization, before QM optimzation and torsion                                                         fitting.

* -opt-b_c_*_angle_postQMOPTprefit.alz          Output from Tinker analyze, before QM optimzation and torsion fitting.

* -sp-b_c_*_angle.com                           Gaussian input for single point calculation.
* -sp-b_c_*_angle.psi4                          Psi4 input for single point calculation.
* -sp-b_c_*_angle.chk                           Unreadable density information from single point calculation.
* -sp-b_c_*_angle.log                           Gaussian or Psi4 output from single point calculation.

* -opt-b_c_*_angle_postQMOPTpostfit.xyz         Tinker XYZ from previous Gaussian log file,after current QM optimzation                                                    and after torsion fitting.

* -opt-b_c_*_angle_postQMOPTpostfit.key         Tinker key using .key_4 parameters but with additional torsion restraints                                                  to fix other torsions during minimization, after QM optimzation and                                                        after torsion fitting.

* -opt-b_c_*_angle_postQMOPTpostfit.xyz_2       Minimized Tinker XYZ, after QM optimzation and after torsion fitting.

* -opt-b_c_*_angle_postQMOPTpostfit.out         Output from Tinker minimization, before QM optimzation and after torsion                                                   fitting.

* -opt-b_c_*_angle_postQMOPTpostfit.alz         Output from Tinker analyze, after QM optimzation and torsion fitting.

* -fit-b-c.png                                    Plot of QM-MM as a function of dihedral angle and the fourier fit.
* -energy-b-c.png                                 Plot of MM1 (prefit),MM2 (postfit),QM and MM1+Fit as a function of                                                         dihedral angle.


### Fragmenter
* datebase.prm                   After all fragments have successfully finished their poltype jobs and then                                    parameters have been transfered to the parent, valence.prms will have the                                     correct format for SMART string (of the fragment) and the corresponding                                       torsion parameters that can be added to the valence.py database.
* structurename.png              2D visualization of input structure with all rotatable bonds with missing                                     torsion parameters shown. Compare this with images of fragments for each                                      highlighted rotatable bond in this image below.
* Below files are in Fragmenter folder
* b_c_Index_i.png                b_c are the middle bond indexes of the torsion a_b_c_d from parent molecule                                   and index i refers to the index in a list of possible combinations of                                         hydrated atoms where we only allow atoms with cut bonds to not be hydrated if                                 the parent atom had a formal charge (not currently checking this when                                         hydrating.If there is only Index_0 this indicated no charged formal atoms                                     from parent.The highlighted rotatable bond corresponds to the b_c highlighted                                 bond in structurename.png (in poltype job top directory).

* FragmentGrowthWithIndex_Bnd_b_c_Index_0.png  Tracks the growth of fragment along while attempting to                                                       minimize differences in Wiberg Bond Order between each fragment                                               and the parent until the tolerance threshold is reached. If                                                   WithIndex in name, this means it will show atoms with their                                                   index numbers. If WithoutIndex in name, it will show atoms                                                    without their index numbers. The last molecule in image is                                                    the molecule that will be used to derive torsion paramaters                                                   for the parent molecule.



* CombintationsWithIndex_Bnd_b_c_Index_0.png   There are different combinations of ways to grow a fragment out                                               and each way is attempted before choosing the minimum                                                         difference in Wiberg bond order between parent bond and                                                       fragment bond. This minimum fragment is added to the                                                          FragmentGrowth image (last fragment on Combinations image).

* parentindextofragindex.txt        Map from parent atom indices to fragment atom indices
* classkeytosmartsposarray.txt      Map from parents torsion key (four type numbers) to the four indices of the                                          torsion atoms in the SMARTS string.
* fragclasskeytoparentclasskey.txt  Map from fragment torsion class key to parent torsion class key
* parentsymclasstofragsymclass.txt  Map from parent type class to fragment type class
* classkeytotorsionindexes.txt      Map from parent torsion key to the torsion atom indices
* classkeytosmarts.txt              Map from parent torsion key to SMART string

### vdw-refinement
* classkeytoatomindexes.txt  Map from parent vdw atom type to atom index locations in SMARTS string
* missingvdw.txt             Identifies which vdw atom indices to be used (missing parameters for)                                                for fragmeneter and vdwfit modules
* -water_a_b.xyz             Dimer XYZ file generated from poltype between water and molecule being                                               parameterized. a is the atom index in the molecule and b is the atom index on                                        the molecule-water dimer.
* -water_a_b_cartesian.xyz   XYZ format of molecule-water dimer. 
* -water_a_b_tinker.xyz      Tinker XYZ format of molecule-water dimer.
* -water_a_b.com             Gaussian input file for geometry optimization.
* -water_a_b.psi4            Psi4 input file for geometry optimization.
* -water_a_b.log             Gaussian or Psi4 output file for geometry optimization.
* -water_a_b.chk             Unreadable electron density from geometry optimization.
* -water_a_b.fchk            Human readable electron density from geometry optimization.
* -water_a_b_X.xyz           Shifted water probe about QM opt minima, where X can be .8,.9,1,1.1,1.2
* -water_a_b_X_sp.com        Gaussian input file for QM SP calculation to obtain counterpoise BSSE corrected energy.
* -water_a_b_X_sp.psi4       Psi4 input file for QM SP calculation to obtain counterpoise BSSE corrected energy.
* -water_a_b_X_sp.log        Gaussian or Psi4 output file for QM SP calculation to obtain counterpoise BSSE                                       corrected energy.
* -water_a_b_X_sp.chk        Unreadable density.
* -water_a_b_X_sp.fchk       Human readable density.
* QM_DATA                    Contains quantum SP energies for each XYZ structure.
* SP.dat                     Contains tinker energies for each XYZ structure.
* INITIAL.PRM                Contains vdw type number being fit, initial radius, intial energy depth, lower bound                                 radius, lower bound energy depth, upperbound radius, upperbound energy depth for fitting
* EnergyVsDistance-water_a_b_Y.png  Energy vs distance plot for sampling the van der waals potential energy surface.                                     Y is the tinker type number being fit.
* QMvsAMOEBA-water_a_b_Y.png        QM vs MM energy plot for sampling the van der waals potential energy surface.                                        Y is the tinker type number being fit.
### Parameter Files
* ModifiedResidue.prm       Modified residue parameters are added to this paramter file. This is originally a copy of 
* amoeba09.prm              Small molecule parameters that databaseparser uses to assign molecules.
* amoebabio18.prm           Contains protein parameters used by molecular dynamics for proteins/DNA.
* amoebabio18_header.prm    Header file added to each key file in poltype.
* externalparameterdatabase.txt External database for adding parameters that databaseparser can assign to molecules.
* mm3.prm                   Parameters from MM3 force field.
* smartstoamoebatypedescrip.txt Dictionary mapping of SMART string to tinker description in amoeba09.prm, used for                                   assigning parameters to molecules in databaseparser.


### Poltype Modules
* symmetry.py                Handles detecting symmetry in input molecule.

* multipole.py               Handles multipole frame definitions and GDMA.
* databaseparser.py          Database for searching valence parameters via SMART string matching to input                                         molecule.
* apicall.py                 Handles external calls for QM job submission on clusters. THE FUNCTION                                               CallExternalAPI IN THIS MODULE NEEDS TO BE MODIFIED FOR SPECIFIC CLUSTER                                             SUBMISSION METHODS.
* optimization.py            Handles Gaussian,Psi4 and Tinker minimization.
* electrostaticpotential.py  Handles electrostatic potential fitting and Gaussian, Psi4 SP calculations.
* torsionfit.py              Handles torsion fitting.
* torsiongenerator.py        Handles torsion parameterization via Gaussian,Psi4,Tinker.
* rings.py                   Handles set up for non-aromatic ring torsion puckering that is used by torsiongenerator                              and torsionfit.
* vdwfit.py                  Handles van der waals paramterization with water dimers.
* fragmenter.py              Handles fragmenting parent molecule, paramterizing torsion and transfering                                           paramters to parent.
* test_poltype.py            Handles unittesting for POLTYPE modules.

### Poltype Examples

