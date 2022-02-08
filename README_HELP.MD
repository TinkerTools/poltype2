## Help

* [Running Poltype](#running-poltype)
* [Input Structure](#input-structure)
* [Most Common Inputs](#most-common-inputs)
* [Geometry Optimization](#geometry-optimization)
* [Distributed Multipole Analysis (DMA) and Electrostatic Single Point (ESP) Inputs](#distributed-multipole-analysis-(dma)-and-electrostatic-single-point-(esp)-inputs)
* [Torsion Parameterization](#torsion-parameterization)
* [Modified Amino Acids](#modified-amino-acids)
* [FAQ](#faq)
* [Tips](#tips)
* [Common Errors](#common-errors)

### Running POLTYPE
* Need two input files to run: *.ini and *.sdf. The ini file sets poltype controls and it specifies the molecule in sdf format to be used.
* Prepare the poltype.ini and sdf file as shown below
* Navigate to directory containing poltype.ini and .sdf file, and run:

```shell
nohup python /path_to_poltype/poltype.py &
```

```final.xyz``` and ```final.key``` are the resulting structure and parameter files you will need.


### Help
```shell
python /path_to_poltype/poltype.py --help
```

### Recommended Protocol (esp for charged molecules)
* For molecules that are neutral and not Zwitterion (total charge can be 0), PSI4 for all QM steps (psi4 is default behaviour).
* For molecules with more than 1 charge center within 3 bond distances away, POLTYPE turns on (implicit solvent) for QM optimization and torsion energy scan (MM/Tinker is set to use Generalized Kirkwood implicit solvation to be consistent)). DMA & ESP fitting of multipoles always in gas-phase. In case PCM is not on automatically for your molecules, you can mannually add options ot turn on PCM (optpcm=True, toroptpcm=True, torsppcm=True). If gaussian is detected in enviorment, then gaussian will be used for PCM computations such as initial geometry optimization, torsion geometry optimization and torsion SP. 
* PSI4 PCM optimzation or single point energy calculation of large molecules may be problematic. In that case, you may consider to use Gaussian (if concentrated charge centers are detected and Gaussian is found in envioronment, then Gaussian is used by default).
* The most expensive calculation in POLTYPE is the ESP (MP2/aug-cc-pvtz). For large molecules, Gaussain requires large memeory and disk, Psi4 is an order og magnitude faster for MP2 (default is psi4).
* Speical cases for zwitter ions
  ** If you need to use Gaussain (PSI4 PCM failed), and do not have suffcient mem/disk, you could lower the basis set to 6-311++G(2d,2p).
  ** another possibility is to (use_gausoptonly. This will use Gaussian for all opt jobs but PSi4 for DMA and ESP QM.  

 
### Input Structure
* Use an SDF file format for the input structure.
* Use totalcharge=charge keyword in .ini file to specify total charge, otherwise poltype will guess based on formal charge from input valence. 
* Bond topologys (input valence) need to be correct for formal charge assignment.
* If you have a radical, then need to add radical to the input SDF file, otherwise be default poltype will set radicals to zero and guess charge based on input valence.

#### Most Common Inputs
```
structure            specifies input structure for poltype, must include total charge information. 
                     Default=None.
numproc              specifies number of processors used for running quantum jobs. 
                     Default=Will be 80% of total available on node
maxmem               specifies max RAM used for running quantum jobs. 
                     Default=Will be 80% of total available on node
maxdisk              specifies max scratch space used for running quantum jobs. 
                     Default=Will be 80% of total available on node

scratchpath          Path to scratch directory
                     Default=/scratch

consumptionratio     Default % of resources to use on node if not specified by user
                     Default=.8

deleteallnonqmfiles  This will delete all non-qm log files (sometimes errors require restarting and deleting the correct files).
                     Default=True

email                Email address to send crash reports to, will send email when poltype job finished.
                     Default=None

debugmode            Cheap QM, just to run through and check for any software bugs.
                     Default=False

jobsatsametime       Parralelize your QM and fragment poltype jobs for local host. 
                     Default=1


externalapi          Specifies path to external API used for submitting quantum jobs on a cluster. Parralelize QM jobs and fragment poltype jobs on cluseter. For fragment jobs, if you wish to make fragment poltype jobs consume less memory then the parent job on the cluster, also add the keyword "jobsatsametime", this will divide parent memory by "jobsatsametime", (number of processors will be rounded to integer).
                     Default=None.

username             Needed for RenLabCluster daemon (externalapi keyword). This will inform the daemon which username to submit jobs for the queuing system.
                     Default=None


fragmentjobslocal Run fragment jobs on localhost even when externalapi keyword is active
                     Default=False

coresperjob          Number of cores for jobs in parralel
                     Default=2

addhydrogens         Adds hydrogens to unprotonated atoms.
                     Default=False

parentjobsatsametime Number of parent molecules at sametime, informs poltype how much resources to divide from available on the host.
                     Default=1

genprotstatesonly    Output protonation states (at pH=7) in mol files(ProtonationState_0.mol, ...) and then exit.
                     Default=False
```

* Note for very large molecules initial opt, DMA and ESP may require large amounts of memory if quantum job consistently crashes

#### Most Common Input Example
__All input arguments are specified in poltype.ini file__
```
structure=methylamine.sdf
numproc=4
maxmem=20GB
maxdisk=100GB
```
* All booleans can switch their truth value, i.e., dontfrag=True or dontfrag=False

####Submitting Quantum Jobs Externally
```
externalapi       specifies path to external API used for submitting quantum jobs on a cluster. 
                  Default=None.

username          Needed for RenLabCluster daemon (externalapi keyword). This will inform the daemon which username to submit jobs for the queuing system.
                  Default=None
```

### Geometry Optimization
```
optbasisset          specifies basis set for geometry optimization. 
                     Default=6-31G*
gausoptcoords        specifies coordinate system used for gaussian geometry optimization. 
                     Default='' 
                     (will default to modified redundant internal coordinates).
freq                 specifies to do a frequency calculation during geometry optimization. Computes true Hessian 
                     and determines if structure is minimum rather than saddle point, but very expensive. 
                     Default=False
optpcm               specifies to use PCM for geometry optimization. 
                     Default=True if there is no net charge, but some atoms have formal charge.
optmaxcycle          specifies the number of max cycles or iterations allowed during optimization.
                     Default=400
optmethod            specifies quantum method for geometry optimziation. 
                     Default=MP2

use_gausgeomoptonly  specifies to only use gaussian for initial geometry optimization and Psi4 for all other QM stepsi
                     Default=False 

dontusepcm           Dont use pcm for any QM calculation.
                     Default=False

use_qmopt_vdw        Use QM opt after MM opt for vdw dimer.
                     Default=False

use_gau_vdw          Use Gaussian for vdw QM if detected in envioronment.
                     Default=False

generateextendedconf Use extended conformation for geometry optimization with torsions restrained.
                     Default=True
```

### Distributed Multipole Analysis (DMA) and Electrostatic Single Point (ESP) Inputs
```
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
```

### Torsion Parameterization
```
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
                  Default=12 if number of torsion parameters to fit is less than 12, else automatically adds more                      points
torspmethod       specifies the SP method used for getting refined multipoles via fitting to quantum potential. 
                  Default=wB97XD
toroptmethod      specifies the opt method used for optimizing the geometry prior to deriving parameters. 
                  Default=wB97XD
toroptbasisset    specifies which basis set to use for geometry optimization during torsion parameterization. 
                  Default=6-31G*
torspbasisset     specifies which basis set to use for SP during torsion parameterization. 
                  Default=6-311++G**
onlyrotbndslist   specifies to only rotate bonds in the list, onlyrotbndlist=1 2, 3 4, 5 6 
                  Default=[] (empty)

tortor            specifies to sample the torsion-torsion potential energy surface for two adjecent and                                rotatable dihedrals
                  Default=False

fitfirsttorsionfoldphase specifies to add a phase fitting parameter to the first cosine term in the torsion fitting                           function
                         Default=False

refinenonaroringtors specifies to pucker and refine torsion parameters on non-aromatic rings.
                     Default=False

maxgrowthcycles      specifies how many growth cycles are allowed by fragmeneter module.
                     Default=4
```

### Fragmenter
```
maxgrowthcycles   specifies the maximum number of growth cycles before stopping fragmenter progtocol.
                  Default=None

dontfrag          specifies to not use fragmenter library for torsion.
                  Default=False

WBOtol            Tolerance for relative differences in Wiberg bond order between parent and fragment. Still                    trying to figure out the "best" tolerance.
                  Default=.001

jobsatsametime    How many jobs (QM or fragment poltype jobs) to run in parralel on single node, will distribute resources evenly to each job.
                  Default=1

fragmentjobslocal   Run fragment jobs on localhost even when externalapi keyword is active
                    Default=False

nonaroringtor1Dscan This keyword will cut a non-aromatic ring and spin the torsion for a 1D scan. Default behaviour is to transfer alkane torsion parameters.
                    Default=False  
```

### Modified Amino Acids
```
modifiedresiduepdbcode    specifies three letter modified amino acid PDB code. 
                          No default.
mutatedsidechain          specifies modified amino acid structure (*.sdf) to replace with existing side chain. 
                          Structure most include side chain of modified AA, connecting to Ca, which is connected
                          to NH2 and COH for protein alignment. 
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
```

### Atom Typing
```
atmidx                    specifies which number to start poltype type numbers from. 
                          Default=401
transferqmequilvalues     specifies to tranfer QM equilbrium angle and bond lengths
                          Default=True
```

### Van der Waals Parameterization.
```
dovdwscan                specifies to not refine initial guess vdw parameters from database assignment.
                         Default=False

onlyvdwatomindex         Only probe one vdw atom index on molecule.

accuratevdwsp            Use CBS extraction on MP2 with counterpoise on basis sets aug-cc-pV[TQ]Z.
                         Default=False

onlyvdwatomlist          List of Van der Waals atom indices to use in VdW refinement.

addlonepairvdwsites      Add dimer ab initio near lone pair sites on target molecule atom.
```

### Misc.
```
keyfiletoaddtodatabase   specifies to add keyfile parameters to amoeba09 parameter file, located in ParameterFiles/.
                         Default=None

deleteallnonqmfiles  This will delete all non-qm log files (sometimes errors require restarting and deleting the correct files).
                     Default=False

dontusepcm           Dont use pcm for any QM calculation.
                     Default=False

email                Email address to send crash reports to.
                     Default=None

inputmoleculefolderpaths   Path to folders containing sdf files for molecules you want to parameterize with poltype. Other settings in poltype.ini file will be written to poltype.ini files in each folder.

debugmode            All QM methods are HF , basis set is MINIX. This is for testing for software crashes, fast QM.
                     Default=False

```


### Minimal Example
Minimal examples are any example folder in Examples/, that contain Symmetry as the first word. Minimal meaning just specifying the structure input file.


### FAQ
* Q)If I want to change QM methods/basis set do I need to restart POLTYPE job from scratch?
  A)No, POLTYPE will look for output files from QM (see README.MANIFEST) to see if it needs to perform a job or continue. So you can delete the outputlog file, change your input options and resubmit POLTYPE job.

### Tips
* For submitting on node clusers, make sure maxdisk and maxmem are set to the lowest amount available out of all of nodes on cluster otherwise risk diskquota / memory error.
* Be careful when changing optimization methods and coordinates without deleting previous XYZ files, DMA and ESP, fchk, cube files, this can cause eletrostaic fitting errors where your grid does not correspond to the potential file from the new coordaintes. 
* If using PSI4 for optimization with high number of optmaxcycles (default is 5 very low), it will be much slower than Gaussian optimization. Additionally even if using lower number of optmaxcycles, if your molecule is chated, by default we use PCM and this will significantly slow down PSI4 optimization. Add keyword "use_gaus" in poltype.ini file.

### Common Errors

* If QM and MM molecular dipole moments are too different, a ValueError will be raised.

* If QM and MM optimized root mean square difference (RMSD) is too high, a ValueError will be raised.
   * If -opt.log coordinates are not the same as .XYZ file coordinates then this is can cause this error.

* If QM and MM root mean square potential difference (RMSPD) is too high, a ValueError will be raised.

* If QM and MM for torsion fitting is too high, a ValueError will be raised. If for some reason you wish to remove an energy pointfrom the fitting, you can locate the corresponding SP QM log file and delete all of the text in it.

* If your Tinker version is less than 8.7, a ValueError will be raised.

* If using modifiedresidues library and a type number for the given atom index and PDB atom label are not found in the library, a ValueError will be raised. 

* If POLTYPE crashes, you do not need to delete any quantum data or tinker data, POLTYPE will recognize that the file exists and attempt to continue with it.

* Error: segmentation violation, this usually means an issue related to memory (not enough RAM or disk). However since we terminate optimization after few cycles (default 5), Gaussian will throw this error in poltype.log. Poltype will ignore this if it also sees 5 opt cycles completed sucessfully. 



### Less Common Errors
*********************************************
* Error! The bond connectivity before and after structure optimization is different. Check to see if hydrogens are other atoms pop off during optimization.
* If your RMSPD is too high and also you have a molecule of high symmetry, please check frame definitions and alter them as needed using POLEDIT and POTENTIAL command in Tinker to change frames/combine multipoles and refit to electrostatic potential. See Symmetry examples under Examples folder. We believe all symmetry cases are already covered by POLTYPE.

* If you are concerned with vdw parameters, you may check again QM dimer (your_molecule-water) structure and energy.


If there are no python errors as shown above, please read *-poltype.log to traceback and see the last command executed. 

* First try and run the command from terminal manually to reproduce the error

* Next have a reference poltype job from the Examples folder that has successfully completed

* Compare the inputs of the command to the inputs in successful poltype job. Look for anything different that might cause the error.

* Try googling issue such as error in Psi4/Gaussian log file

### Modified Amino Acid Errors:
* Typically errors here that appear as undefined types or missing parameters in the final key files are a result of inconsistency between the atom label and residue label in the input PDB file and the atom label, residue label in the library lib.bio18_conv1.txt used to convert atomlabel, residue label to type numbers.

* Delete any atoms at end of PDB not a part of a full residue (like an alpha carbon at the very end).

* If the number of expected atoms connected to an atom type are not found, it is possible either the incorrect PDB atom label was used for that protonation state or there is a missing/extra hydrogen.

* If restarting poltype job after error delete the ModifiedPrm.prm file in poltype directory before restarting.

### ForceBalance Poltype Wrapper
 
* Install modified forcebalance for Tinker9 GPU dynamics
```
cd ~
git clone https://github.com/bdw2292/forcebalance.git
```

* Make sure to copy files from modified forcebalance to conda site packages (for example mine are located /home/bdw2292/miniconda3/envs/FBTest/lib/python3.8/site-packages/forcebalance/) 
```
cp ~/forcebalance/src/*py /home/bdw2292/miniconda3/envs/amoebamdpoltype/lib/python3.6/site-packages/forcebalance/
cp ~/forcebalance/src/data/*py /home/bdw2292/miniconda3/envs/amoebamdpoltype/lib/python3.6/site-packages/forcebalance/data/

```

* Example input
```
poltypepathlist=/home/bdw2292/PoltypeJobs/HalogenPaperMolecules/CH3Cl,/home/bdw2292/PoltypeJobs/HalogenPaperMolecules/CH2Cl2,/home/bdw2292/PoltypeJobs/HalogenPaperMolecules/CHCl3 
csvexpdatafile=ExpHalogenData.csv
vdwtypeslist=401 402 403, 401 402 403 , 401 402 403  
liquid_prod_time=.5                 
gas_prod_time=.5
WQ_PORT=9013
```
* Example CSV File in ForceBalanceWrapperExamples folder
* Only use WQ_PORT for parralel jobs, otherwise default is serial on single node
* Dont let property values be blank for some TP points (need value for FB to read).
* If using csvexpdatafile, please use correct units as below and the column headers in example csv file, make sure the molecule name under the Name column are same as poltype molecule folder name
* By default all types in key files are reassigned new type numbers so there is no overlap while fitting together. So input types will be converted internally.
* temperature_list, pressure_list and density_list are required.
* Density units are Kg/m^3
* Enthalpy of vaporization units kJ/mol
* Thermal expansion units 1/K
* Isothermal compressibility units are 1/bar
* Isobaric heat capacity units are kJ/mol.K
* pressure units are atm
* Experimental data can be found at https://webbook.nist.gov/ or https://www.aiche.org/dippr
* QM energies from dimers in poltype (normalized) >=15 kcal/mol are removed and not use for forcebalance.
```
poltypepathlist                                      Comma seperated list of Paths to completed poltypejobs
vdwtypes                                             Comma seperated list of vdw types lists to parameterize
WQ_PORT                                              Port number for work queue to connect to workers on nodes
liquid_equ_steps                                     Equilibrium steps MD 
liquid_prod_steps                                    Production dynamic steps MD
liquid_timestep                                      Time step for liquid MD
liquid_interval                                      How often frames written out, in ps
gas_equ_steps                                        Equilibrium steps MD 
gas_prod_steps                                       Production dynamic steps MD
gas_timestep                                          
gas_interval                                         How often frames written out, in ps
md_threads                                           CPU cores
liquid_prod_time                                     Total time liquid sim
gas_prod_time                                        Total time gas sim
debugmode                                            Sets MD to very short time for testing ForceBalance (.001 ns)
```
