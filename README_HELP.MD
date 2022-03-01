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
* The most expensive calculation in POLTYPE is the ESP (MP2/aug-cc-pvtz). For large molecules, Gaussain requires large memeory and disk, Psi4 is an order of magnitude faster for MP2 (default is psi4).
* PCM is turned on by default when zwitterions are detected
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


fragmentjobslocal    Run fragment jobs on localhost even when externalapi keyword is active
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

#### Minimal Input Example
__All input arguments are specified in poltype.ini file__
```
structure=methylamine.sdf
```
* All booleans can switch their truth value, i.e., dontfrag=True or dontfrag=False

#### Submitting Quantum Jobs Externally
```
externalapi       specifies path to external API used for submitting quantum jobs on a cluster. 
                  Default=None.

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


numespconfs       specifies number of conformations to be used for ESP potential fitting. Rdkit conformors are evaluated via MMFF94 energies and numespconfs-1 most favorable conformations are picked to add to fitting (with first being most extended conformation). 
                  Default=1
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


### Atom Typing
```
atmidx                    specifies which number to start poltype type numbers from. 
                          Default=401
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
sleeptime                Number of seconds to wait before each loop iteration in job submission queue.
                         Default=5

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

* If POLTYPE crashes, you do not need to delete any quantum data or tinker data, POLTYPE will recognize that the file exists and attempt to continue with it.

* Error: segmentation violation, this usually means an issue related to memory (not enough RAM or disk).  


### ForceBalance Poltype Wrapper
 

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
vdwprmtypestofit                                     Sets allowed parameter types for all vdw parameters to fit for vdw (radius S, depth  T, shift D), default=S,T,D . Default makes all vdw parameter types able to be fit.
vdwtypestoeval                           For parameters being fixed, this tells poltype which ones to fix (while the first instance of parameters is fit), default all parameters are evaluated besides for instance vdwtypestoeval=S,T,D. 


```

* View ForceBalance results
```
conda activate pymolenv
python /home/bdw2292/poltype2/PoltypeModules/plotFBresults.py /home/bdw2292/FBJobs/
 
```
