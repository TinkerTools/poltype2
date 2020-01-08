HELP
*********************************************
python /path_to_poltype/poltype.py --help

Unittesting
*********************************************
python /path_to_poltype/poltype.py --unittest

Running POLTYPE
*Prepare the poltype.ini file as shown below first*
*Navigate to directory containing poltype.ini and .sdf file*
*********************************************
python /path_to_poltype/poltype.py

Input Structure - SDF File
*********************************************
After the bond block, you can specify formal charge. Just put the total charge on one atom.

The following example means there is 1 atom with a formal charge: atom 2 has a +1 charge

M  CHG  1   2   1

Most Common Inputs
********************************************* 
structure         specifies input structure for poltype, must include total charge information. 
                  Default=None.
numproc           specifies number of processors used for running quantum jobs. 
                  Default=1.
maxmem            specifies max RAM used for running quantum jobs. 
                  Default=700MB 
maxdisk           specifies max scratch space used for running quantum jobs. 
                  Default=100GB

Most Common Input Example
**All input arguments are specified in poltype.ini file**
*********************************************
structure=methylamine.sdf
numproc=4
maxmem=20GB
maxdisk=100GB


Submitting Quantum Jobs Externally
********************************************* 
bashrcpath        specifies path for bashrc to source for running quantum jobs, used when calling external API 
                  for submitting quantum jobs. 
                  Default None.
externalapi       specifies path to external API used for submitting quantum jobs on a cluster. 
                  Default=None.

Geometry Optimization
********************************************* 
optbasisset       specifies basis set for geometry optimization. 
                  Default=6-31G*
gausoptcoords     specifies coordinate system used for gaussian geometry optimization. 
                  Default='' 
                  (will default to modified redundant internal coordinates).
freq              specifies to do a frequency calculation during geometry optimization. Computes true Hessian 
                  and determines if structure is minimum rather than saddle point, but very expensive. 
                  Default=False
optpcm            specifies to use PCM for geometry optimization. 
                  Default=True if there is a net charge.
optmaxcycle       specifies the number of max cycles or iterations allowed during optimization.
                  Default=400
optmethod         specifies quantum method for geometry optimziation. 
                  Default=MP2

Distributed Multipole Analysis (DMA) and Electrostatic Single Point (ESP) Inputs
********************************************* 
use_gaus          specifies to use Gaussian for quantum geometry optimzation and SP calculations. 
                  Default=False
use_gausoptonly   specifies to use Gaussian for Geometry optimization only and Psi4 for SP calculations. 
                  Default=False
qmonly            specifies to run quantum for geometry optimization, DMA SP and ESP SP then quit. 
                  Default=False.
dmamethod         specifies which quantum method to use for distributed multipole analysis (DMA). 
                  Default=MP2
espmethod         specifies which quantum method to use for refined multipoles. 
                  Default=MP2
dmabasisset       specifies basis set to use for quantum that generates initial multipoles. 
                  Default=6-311G**
espbasisset       specifies basis set to use for quantum intended for refined multipoles. 
                  Default=6-311++G(2d,2p)

Torsion Parameterization
*********************************************       
rotalltors        specifies poltype to spin all rotatable bonds instead of transferring torsion parameters if 
                  found in database. 
                  Default=False
torsppcm          specifies poltype to use PCM for torsion SP. 
                  Default=True if there is net charge
toroptpcm         specifies poltype to use PCM for torsion optimization. 
                  Default=True if there is net charge.
dontdotor         specifies poltype to skip torsion spinning. 
                  Default=False
dontdotorfit      specifies poltype to skip torsion fitting. 
                  Default=False
torsionrestraint  specifies how many kcals/mol to restrain other torsions around rotatable bond to be at 
                  current dihedral angle. 
                  Default=.1
foldnum           specifies the number of fourier terms used in torsion fitting.
                  Default=3
tordatapointsnum  specifies the number of quantum jobs to sample along the dihedral angle for all torsions.
                  Default=12 if number of torsion parameters to fit is less than 12, else automatically adds more points
fitrotbndslist    specifies poltype to derive torsion parameters for all bonds
                  specified as fitrotbndslist=1 2, 3 4, 5 6
                  Default=[] (empty)
torspmethod       specifies the SP method used for getting refined multipoles via fitting to quantum potential. 
                  Default=MP2
toroptmethod      specifies the opt method used for optimizing the geometry prior to deriving parameters. 
                  Default=HF
toroptbasisset    specifies which basis set to use for geometry optimization during torsion parameterization. 
                  Default=6-31G*
torspbasisset     specifies which basis set to use for SP during torsion parameterization. 
                  Default=6-311++G**
onlyrotbnd        specifies to only rotate single bond, onlyrotbnd=1,2 . 
                  Default=[] (empty)

Modified Amino Acids
*********************************************
modifiedresiduepdbcode    specifies three letter modified amino acid PDB code. 
                          No default.
mutatedsidechain          specifies modified amino acid structure (*.sdf) to replace with existing side chain. 
                          Structure most include immediate backbone for protein alignment. 
                          No default.
mutatedresiduenumber      specifies which amino acid residue will be replaced with mutated side chain. 
                          No default.
unmodifiedproteinpdbname  specifies the PDB (with protonation states already assigned), with un-modified 
                          amino acids. 
                          No default.
modifiedproteinpdbname    specifies the PDB with modified amino acids, if not replacing existing side chain, or 
                          this was already done. 
                          No default.
amoebabioprmpath          speicfies the path to amoebabio18.prm file for stitching parameters from existing 
                          protein/NA parameters to mutated parameters. 
                          No default.

Atom Typing
********************************************* 
atmidx   specifies which number to start poltype type numbers from. 
         Default=401
uniqidx  specifies to use atom indexes as atom type numbers rather than averaging types
         Default=False


*********************************************
Check your results!
*********************************************
ttt.xyz and ttt.key are the resulting structure and parameter files you will need.

If QM and MM molecular dipole moments are too different, a ValueError will be raised.

If QM and MM optimized root mean square difference (RMSD) is too high, a ValueError will be raised.

If QM and MM root mean square potential difference (RMSPD) is too high, a ValueError will be raised.

If QM and MM for torsion fitting is too high, a ValueError will be raised.

For molecules of high symmetry, please check frame definitions and alter them as needed using POLEDIT and POTENTIAL command in Tinker to change frames/combine multipoles and refit to electrostatic potential. See Symmetry examples under Examples folder.

If you are concerned with vdw parameters, you may check again QM dimer (your_molecule-water) structure and energy.


