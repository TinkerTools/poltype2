## Advanced Usage

### Advanced Program Usage Parameterization

[Required Inputs Parameterization](#required-inputs-parameterization)

[Controlling Resource Consumption](#controlling-resource-consumption)

[Protonation State](#protonation-state)

[Controlling Quantum Packages](#controlling-quantum-packages)

[Atom Typing](#atom-typing)

[Geometry Optimization](#geometry-optimization)

[Distributed Multipole Analysis and Electrostatic Single Point](#distributed-multipole-analysis-and-electrostatic-single-point)

[Database Matching](#database-matching)

[Torsion Parameterization](#torsion-parameterization)

[Fragmenter](#fragmenter)

[Van der Waals Parameterization](#van-der-waals-parameterization)

### ForceBalance Poltype Wrapper

[ForceBalance Poltype Wrapper Example](#forcebalance-poltype-wrapper-example)

[ForceBalance Poltype Wrapper Keywords](#forcebalance-poltype-wrapper-keywords)

[ForceBalance Poltype Wrapper Results](#forcebalance-poltype-wrapper-results)


### Advanced Program Usage Molecular Modelling

[Protein Preparation Inputs](#protein-preparation-inputs)

[Docking Inputs](#docking-inputs)

### Advanced Program Usage Molecular Dynamics

[Required Inputs Molecular Dynamics](#required-inputs-molecular-dynamics)

[Binding Free Energy](#binding-free-energy)

[Solvation Free Energy](#solvation-free-energy)

[Neat Liquid Simulation](#neat-liquid-simulation)

[Generic Molecular Dynamics Inputs](#generic-molecular-dynamics-inputs)

[Available Integrators Barostats and Thermostats](#available-integrators-barostats-and-thermostats)

[Box Setup](#box-setup)

[Minimization](#minimization)

[Restraints](#restraints)

[Equilibriation](#equilibriation)

[Production Dynamics](#production-dynamics)

[Free Energy](#free-energy)

[Mutation](#mutation)


#### Notes
* All booleans can switch their truth value, i.e., addhydrogens=True or addhydrogens=False
* For very large molecules initial opt, DMA and ESP may require large amounts of memory if quantum job consistently crashes
* Defaults are chosen specifically for maximizing accuracy and efficiency (without sacrificing accuracy)
* Be careful when changing optimization methods and coordinates without deleting previous XYZ files, DMA and ESP, fchk, cube files, this can cause eletrostaic fitting errors where your grid does not correspond to the potential file from the new coordaintes. 



#### Required Inputs Parameterization
``structure`` specifies input structure for poltype as an SDF file. Must include bond orders consistent with desired total charge. Default is None.



#### Controlling Resource Consumption
``numproc`` specifies the number of processors allowed. Default will be 80% of total available on node.

``maxmem`` specifies the maximum amount of RAM allowed. Default will be 80% of total available on node.


``maxdisk`` specifies maximum scratch space allowed. Default will be 80% of total available on node.


``consumptionratio`` specfies ratio of resources to use on node if not specified by user. Default is .8.

``coresperjob`` specifies number of cores per job for jobs in parralel. Default is 2.

``jobsatsametime``  specifies how many fragment jobs are allowed to run in parralel on node. Default is int(numproc/coresperjob) for input poltype job but 1 for fragment job.

``coresperjob`` specifies number of cores per job for jobs in parralel. Default is 2.

``parentjobsatsametime`` specifies number of input poltype jobs running on same host, rhis informs poltype how much resources to divide from available on the host. Default is 1.

#### Protonation State
``addhydrogens`` specifies to adds hydrogens to unprotonated atoms. Default is False.


``genprotstatesonly``    specifies to output protonation states (at pH=7) in mol files(ProtonationState_0.mol, ...) and then exit. Default is False.


#### Controlling Quantum Packages
``use_gaus`` specifies to use Gaussian for all QM calculations. Default is False.


#### Atom Typing

``atmidx``  specifies which number to start poltype type numbers from. Default is 401.


``usesymtypes`` specifies to use symmetry typing to reduce and average atom types. Default is True.


``indextotypefile`` specfies to use a file that will read in custom atom index to type number map. Default is None. Example ``indextotypefile=indextotypefile.txt``. For each line in the file, the first item is atom index (starting from 1) and second item is the desired type number.

#### Geometry Optimization

``optbasisset`` specifies basis set for geometry optimization. Default is 6-31G*.

``optmethod`` specifies quantum method for geometry optimziation. Default is MP2.
   

``optpcm`` specifies to use PCM for geometry optimization. Default is True if zwitterion is detected, False if not.

``use_gausoptonly``  specifies to only use gaussian for initial geometry optimization and Psi4 for all other QM steps. Default is False. 

``generateextendedconf`` specifies to use extended conformation for geometry optimization with torsions restrained. Default is True.


#### Distributed Multipole Analysis and Electrostatic Single Point


``dmamethod``  specifies which quantum method to use for distributed multipole analysis (DMA). Default is MP2.
                  
``dmabasisset`` specifies basis set to use for quantum that generates initial multipoles. Default is 6-311G**.         
                  
``espmethod``   specifies which quantum method to use for electrostatic fitting. Default is MP2.

``espbasisset`` specifies basis set to use for electrostatic fitting. Default is aug-cc-pvtz.

``numespconfs`` specifies number of conformations to be used for ESP potential fitting. Rdkit conformors are evaluated via MMFF94 energies and numespconfs-1 most favorable conformations are picked to add to fitting (with first being most extended conformation). Default is 1.

``espextraconflist`` specifies list of conformors to add to extended conformation for ESP fitting. Default is to just use extended conformation for ESP fitting. Files can be SDF or any format as long as atom order is same as input molecule atom order.

``indextompoleframefile`` specifies file of multipole frame definitions for each line (exactly how frames are defined in first key file). This will be used instead of poltype frame detection if specified.


#### Database Matching

``writeoutmultipole`` specifies to write out multipole terms in pre-fit vdw key (from database matches). Default is True.

``writeoutbond`` specifies to write out bond terms in pre-fit vdw key (from database matches). Default is True.

``writeoutangle`` specifies to write out angle terms in pre-fit vdw key (from database matches). Default is True.

``writeoutstrbnd`` specifies to write out strbnd terms in pre-fit vdw key (from database matches). Default is True.

``writeoutopbend`` specifies to write out opbend terms in pre-fit vdw key (from database matches). Default is True.

``writeoutvdw`` specifies to write out vdw terms in pre-fit vdw key (from database matches). Default is True.

``writeoutpolarize`` specifies to write out polarize terms in pre-fit vdw key (from database matches). Default is True.

``writeouttorsion`` specifies to write out torsion terms in pre-fit vdw key (from database matches). Default is True.

``inputkeyfile`` specifies a keyfile to append parameters to pre-fit vdw key (from database matches). Default is None.

#### Torsion Parameterization

``rotalltors``  specifies poltype to spin all rotatable bonds instead of transferring torsion parameters if found in database. Default is False.
                  
``torsppcm``    specifies poltype to use PCM for torsion SP. Default is True if zwitterion is detecter, False otherwise


``toroptpcm``   specifies poltype to use PCM for torsion optimization. Default is True if zwitterion is detecter, False otherwise.
                  
``dontdotor``   specifies poltype to skip torsion parameterization. Default is False.
                  

``torsionrestraint``  specifies how many kcals/mol/rad^2 to restrain other torsions around rotatable bond to be at current dihedral angle. Default is 1641.4 (.5 kcal/mol/deg^2).

``foldnum``           specifies the number of fourier terms used in torsion fitting. Default is 3.

``tordatapointsnum``  specifies the number of quantum jobs to sample along the dihedral angle for all torsions. Default is 12 if number of torsion parameters to fit is less than 12, else automatically adds more points.

``toroptmethod`` specifies the opt method used for torsion parameterization. Default is wB97XD.


``toroptbasisset``    specifies which basis set to use for geometry optimization during torsion parameterization. Default is 6-31G*.

``torspmethod`` specifies the SP method for torsion parameterization. Default is wB97XD.

``torspbasisset``  specifies which basis set to use for SP during torsion parameterization. Default is 6-311+G*.
                  

``onlyrotbndslist``   specifies to only rotate bonds in the list, onlyrotbndlist=1 2, 3 4, 5 6 . Default has no value.


``dontrotbndslist``   specifies to not rotate bonds in the list, onlyrotbndlist=1 2, 3 4, 5 6 . Useful for modified AA lib, to tell poltype not to rotate backbone dihedrals. Default has no value.


``tortor``            specifies to sample the torsion-torsion potential energy surface for two adjecent and rotatable dihedrals. Default is False.

``onlyrottortorlist`` specifies list of tortors to rotate. "onlyrottortorlist=b1 c1 d1 , b2 c2 d2" where b1,c1,d1 belong tor tortor a1,b1,c1,d1,e1. Default has no value.
  

``refinenonaroringtors`` specifies to pucker and refine torsion parameters on non-aromatic rings. Default is False.


``nonaroringtor1Dscan`` specifies to cut a non-aromatic ring and spin the torsion for a 1D scan. Default behaviour is to transfer alkane torsion parameters.
                    
``fitfirsttorsionfoldphase`` specifies to add a phase fitting parameter to the first cosine term in the torsion fitting                           function. Default is False.

#### Fragmenter

``maxgrowthcycles``   specifies the maximum number of growth cycles before stopping fragmenter progtocol. Default is 4.

``dontfrag``          specifies to not fragment input molecule. Default is False.

``WBOtol``            specifies the tolerance for relative differences in Wiberg bond order between parent and fragment. Default is .05.



#### Van der Waals Parameterization

``dovdwscan``                specifies to not refine initial guess vdw parameters from database assignment with ab initio dimer data. Default is False.

``onlyvdwatomindex`` specifies to probe one vdw atom index on molecule. Default has no value.

``accuratevdwsp``   specifies to use CBS extraction on MP2 with counterpoise on basis sets aug-cc-pV[TQ]Z. Default is False.

``onlyvdwatomlist`` specifies list of Van der Waals atom indices to use in VdW refinement. Default has no value.

``fixvdwtyperadii`` specifies list of vdwtypes to fix radius in fitting (can let other vdw parameters fit). Default has no value.




### ForceBalance Poltype Wrapper
 
 
#### ForceBalance Poltype Wrapper Example
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
* Dont let property values be blank for some TP points (need value for FB to read).
* If using ``csvexpdatafile``, please use correct units as below and the column headers in example csv file, make sure the molecule name under the Name column are same as poltype molecule folder name
* Density units are Kg/m^3
* Enthalpy of vaporization units kJ/mol
* pressure units are atm
* QM energies from dimers in poltype (normalized) >=10 kcal/mol are removed and not use for forcebalance.

#### ForceBalance Poltype Wrapper Keywords


``poltypepathlist`` specifies comma seperated list of Paths to completed poltypejobs.

``vdwtypes``        specifies comma seperated list of vdw types lists to parameterize.

``WQ_PORT``         specifies port number for work queue to connect to workers on nodes. Only used for job parralelization.

``liquid_equ_steps`` specifies steps for liquid equilbriation dynamics.

``liquid_prod_steps`` specifies steps for liquid production dynamics.

``liquid_timestep``  specifies time step for liquid dynamics.

``liquid_interval`` specifies how often frames written out, in ps for liquid dynamics.

``gas_equ_steps``   specifies steps for gas equilbriation dynamics.

``gas_prod_steps``  specifies steps for gas production dynamics.

``gas_timestep``   specifies time step for gas dynamics.  

``gas_interval``   specifies how often frames written out, in ps for gas dynamics.

``md_threads``  specifies number of CPU threads to use.

``liquid_prod_time`` specifies total time for liquid production dynamics (ns).

``gas_prod_time``   specifies total time for gas production dynamics (ns).

``liquid_equ_time`` specifies total time for liquid equilibriation dynamics (ns).

``gas_equ_time``    specifies total time for gas equilibriation dynamics (ns).


``vdwprmtypestofit`` specifies allowed parameter types for all vdw parameters to fit for vdw (radius S, depth  T, shift D), default=S,T,D . Default makes all vdw parameter types able to be fit.

``vdwtypestoeval``  specifies which vdw parmaeters to fix (while the first instance of parameters is fit). Default all parameters are evaluated besides for instance vdwtypestoeval=S,T,D. 

``qmrelativeweight`` specifies relative weight between QM and Liquid target. Default is .05. 

``liqrelativeweight`` specifies relative weight between Liquid and QM target. Default is 2.5. 

``enthalpyrelativeweight`` specifies relative weight between enthalpy and other liquid properties. Default is 1. 

``densityrelativeweight`` specifies relative weight between density and other liquid properties. Default is .01.

``targetdensityerror`` specifies target density error in kg/m^3 for ForceBalance fitting. Default is 10 kg/m^3.

``targetenthalpyerror`` specifies target enthalpy error in kcal/mol for ForceBalance fitting. Default is .24 kcal/mol .  

``relaxFBbox`` specifies to equilbriate neat liquid box at room temperature for ForceBalance input files (.5 ns and time step of 2 fs). Default is False.

#### ForceBalance Poltype Wrapper Results

``forcebalancejobsdir=Path_to_force_balance_jobs`` specifies to grab ForceBalance output files, parse and generate plots and tables.


### Advanced Program Usage Molecular Dynamics


#### Protein Preparation Inputs

``pdbcode`` specifies to take PDB code, downloads the crystal from RCSB and used Modeller to find missing residues and model them. Fast optimizaiton options are used by default. Missing residue information is found in REMARKS section of PDB which poltype reads. Default pdbcode is None.


``usepdb2pqr`` specifies to use PDB2PQR program to compute pKa values for titratable residues via propka then will protonate the PDB for you and quit program. Default is False.

#### Docking Inputs

``dockgridcenter`` specifies center in protein pocket to place ligand for docking. Default will use center of ligand.

``dockgridsize`` specifies how far grid extends for docking. Default is dockgridsize=20,20,20 .

``vinaexhaustiveness`` specifies how exaustive to search on grid for AutoDock Vina and vinardo. Default is 32.

``nposes`` specifies number of poses for docking. Default is 10.

``usevina`` specifies to use AutoDock Vina. Default is False.

``usead4`` specifies to use AutoDock 4. Default is False.

``usevinardo`` specifies to use Vinardo scoring function. Default is False.

``usegold`` specifies to use GOLD. Default is False. 

``gridspacing`` specifies grid density for AutoDock 4. Default is .4 angstroms.


``ecrexpect`` specifies expected value of exponential distribution for exponential conensus ranking. This roughly controls how many molecules will be taken into consideration for consensus. Default is 10.


#### Required Inputs Molecular Dynamics
``ligandxyzfilenamelist``  specifies the list of input tinker XYZ filename for the ligands. Default has no value.

``keyfilenamelist`` specifies the list of input tinker key file for the ligand parameters. Default is None.


#### Binding Free Energy

``binding`` specifies to compute solvation and complexation free energies and combine for binding free energy. Default is False.


``complexedproteinpdbname``     specifies complexed protein PDB file, needed for extracting the ligand and converting to tinker XYZ file. Make sure no atoms are missing!!. The atom order between ligand in this file and ligandxyzfile (tinker XYZ) needs to be same. Default has no value.


``receptorligandxyzfilename``   specifies tinker XYZ file contains ligand and protein. This is automatically generated when the above complexedproteinpdbname is given. Default has no value.
                                      


                                      
#### Solvation Free Energy
``solvation``   specifies to do solvation simulations. Default is False.

``salthfe``     specifies whether to compute salt hydration free energy.  Default is False.

#### Neat Liquid Simulation

``neatliquidsim``  specifies to do neat liquid simualtion. Default is False.

#### Generic Molecular Dynamics Inputs 

``annihilateligandxyzfilenamelist``  specifies the list of input tinker XYZ filename for the ligands that will be disappeared. Default is same as ``ligandxyzfilenamelist``.


``integrator``  specifies the MD integrator method. Default is RESPA.

``thermostat``  specifies the thermostat method during NVT MD. Default is BUSSI.

``barostatmethod``  specifies the barostat method used during NPT. Default is montecarlo.

``vdwcutoff``   specifies the cutoff radius for van der Waals interactions. Units in angstroms. Default is 12.

``ewaldcutoff`` specifies the cutoff radius for electrostatic interactions. Units in angstroms. Default is 7.

``polareps``    specifies the Grid spacing for electrostatic potential surface. Default is .00001.

``templateligandxyzfilename``   specifies input template tinker XYZ filename for aligning ligandxyzfilenamelist (first item in list) to templateligandxyzfilename. Program will quit after alignment. Default has no value. 

``externalapi`` specifies external script path for submitting MD jobs. Used with submitjobs.py. If user creates their own submission script, must edit submitjobs.py to work with it. Default is to submit jobs on local node. The submission script needs to be able to handle
correct assignment of tinker envioronments on whichever CPU/GPU cluster it operates on. If using RenLabCluster please read the file on https://github.com/bdw2292/RenLabCluster. Default has no value.


``numproc`` specifies the number of threads to use for CPU dynamics. Default is 1.

``prmfilepath`` specifies the path to parameters. Default is amoebabio18.prm.

``simulationstostopfolderpath``  specifies a directory containing subdirectories of input jobs, recursive search will and *.end files to each folder with an *.arc file to stop simulations. Upon restarting program all *.end files that are in folders are automatically removed. Default has no value.

#### Available Integrators Barostats and Thermostats

##### Integrators

* Verlet                    ``integrator=verlet``
* Beeman                    ``integrator=beeman``
* RESPA                     ``integrator=respa``
* Stochastic                  ``integrator=stochastic`` 
* Nose-Hoover NPT             ``integrator=nose-hoover``


##### Barostats

* Berendsen                ``barostat=berendsen``
* Monte-Carlo              ``barostat=montecarlo``
* Lagevin                  ``barostat=Langevin``

##### Thermostats

* Bussi-Parrinello          ``thermostat=bussi``
* Berendsen                 ``thermostat=berendsen``
* Andersen Stochastic       ``thermostat=andersen``
* Nose-Hoover               ``thermostat=nose-hoover``






#### Box Setup

``density`` specifies the density of neat liquid box in kg/m^3. Default has no value.


``aaxis``   specifies the fixed box size for aaxis of PBC box. If not used, will be determined automatically. Default has no value.

                                        
``baxis``   specifies the fixed box size for baxis of PBC box. If not used, will be determined automatically.  Default has no value.
                                        

``caxis``   specifies the fixed box size for caxis of PBC box. If not used, will be determined automatically. Default has no value.   
                                        

``listofsaltcons``  specifies the salt concentration for physiological systems. Units in mM. Default is [KCl]=100.


``usepreequilibriatedbox``  specifies to use preequilibriated box and then cut default equilibriation time in half. Default is False.

``boxonly`` specifies to stop program after making box. Default is False.           


``addphysioions`` specifies to add physiological ions to box system to mimic envioronment inside cell. Default is KCl at 100mM concentration. 

#### Minimization

``restrainpositionconstant`` specifies the restraint position force constant during minimization for freezing the ligand and protein. Units in kcal/mol/Ang^2.  
                                        
``loosemincriteria``    specifies the convergence criteria for minimization, loose criteria. Units kcal/mol. Default is 10.

``minonly`` specifies to quit after minimization is complete. Default is False.

``minfinished`` specifies to skip minimization (for debugging). Default is False. 

#### Restraints

``restrainatomgroup1``  specifies the group of atoms (protein or ligand) used for tinker group restraints. If not specified, will be automatically determine (for group restraints during equilbriation). Default has no value.


``restrainatomgroup2``  specifies the group of atoms (protein or ligand) used for tinker group restraints. If not specified, will be automatically determine (for group restraints during equilbriation). Default has no value.

``flatbotrest``     specifies whether to use flatbottom potential with group restraints. Default is True.


``distancerestraintconstant``  specifies the restraint constant for restrain-distance keyword during production dynamics. Units in kcal/mol/Ang^2. Default is 10.
                         


#### Equilibriation

``equilrestrainsphereradius``   specifies the sphere radius for protein atoms to restrain position during equilibriation. Units angstroms. Default is 2.

``equiltimeNVT`` specifies the total time in NVT ensemble during equilbriation (for all temperatures combined except last temperature which is repeated by default). Units in ns. Default is 4.


``equilrestlambdascheme``   specifies the group restraint force constants used during equilibriation (NVT). Default is 5,2,1,.5,.1,0.


``equilibriatescheme``  specifies the temperature used for each NVT step during equilibriation. Units in K. Default is 50,100,150,200,300,300

``equilwritefreq``  specifies how often to write out new frame during equilbriation MD. Units in ps. Default is 100.

``equiltimestep``   specifies the MD time step during equilbriation. Units in fs. Default is 2.


``equiltimeNPT``    specifies the total time in NPT ensemble (will determine average box length if using NVT production dynamics). Units in ns. Default is 1.


``equiltimeionNPT``    specifies the total time in NPT ensemble for ion box when computing salt hydration free energy. Units in ns. Default is 1.


``lastNVTequiltime`` specifies the last NVT equilibriation time longer than the previous ones. The previous ones are computed as equiltimeNVT/number_temperatures. Default is .5.


``equilonly``   specifies to stop program after equilibriation and setup for production dynamics. Default is False.


``equilfinished`` specifies to skip equilibration (for debugging). Default is False.

#### Production Dynamics

``productiondynamicsNVT``   specifies whether or not to use NVT ensemble during production dynamics. If on then will use NPT after equilbriating to determine average box size. Default is False.

                                    
``restlambdascheme``    specifies the restraint lambda scheme during production dynamics. Default is 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.

``vdwlambdascheme``     specifeis the Van der Waals lambda scheme during production dynamics. Default is 0,.45,.52,.56,.58,.6,.62,.64,.67,.7,.75,.8,.85,.9,.95,1,1,1,1,1,1,1,1,1,1,1,1.

``estatlambdascheme``   specifies the electrostatics lambda scheme during production dynamics. Default is 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1.

``proddynwritefreq``    specifies how often to write out new frame during production dynamics MD. Units in ps. Default is 2.


``proddyntime``  specifies the total time to run production dynamics (in NPT by default) for each lambda. Units in ns. Default is 5 for solvation and 10 ns for binding computations.

``proddyntimestep`` specifies the MD time step during production dynamics. Units in fs. Default is 2.

``checktraj``   specifies to check each trajectory for any corruption and delete if found. Default is False (slow to read each arc file).

``prodmdfinished`` specifies to skip production dynamics (for debugging). Default is False.

#### Free Energy

``barinterval``  specifies how many times to recompute BAR over length of full production trajectory with evenly spaced chunks over the trajectory. Default is 5.

``perturbkeyfilelist``  specifies a comma seperated list of keyfiles to add to end state and perturb with BAR. Right most item in the list is the end state. By default, arc file from nearestneighboring trajectory is used if fep is turned on, otherwise running MD for new state is default. Default has no value.

``numinterpolsforperturbkeyfiles``  specifies a comma seperated list of number of interpolations per key file in perturbkeyfilelist (they should be same length). Use 1 for keyfiles that need only one interpolation (just between neighboring states), otherwise parameters will be interpolated via step=(b-a)/n where n is number of interpolations for that keyfile. Default is 1 for all keys in perturbkeyfilelist.

``expfreeenergy``  specifies the experimental free energy for comparison in output tables. Default has no value.


``fep`` specifies to copy nearest MD trajectory instead of running MD for each new keyfile to be perturbed. Used with perturbkeyfilelist keyword. Default is False.


``pathtosims``  specifies path to folder containing subfolders of simulations for combining simulation data into one spreadsheet. Will exit program after combining data. Default has no value. 

``extractinterforbinding`` specifies to extract intermoleculer contribution to free energy for each BAR step. By default, annihilation is used which causes inter and intra molecular components to be mixed for binding free energy (gas phase is not computed since it would cancel in thermodynamic cycle). This keyword will extract the ligand dynamics and generate an ARC file for the ligand, run BAR and compute the exact intermolecular component to the free energy and add to spreadsheet.

``barfilesfinished`` specifies to skip BAR (for debugging). Default is False.

#### Mutation

``bgnstatekey``  specifies the first state key file for mutation. Default has no value.

``endstatekey``  specifies the end state key file for mutation. Default has no value.

``bgnstatexyz`` specifies the first state xyz file for mutation. Default has no value.

``endstatexyz`` specifies the end state xyz file for mutation. Default has no value.

``bgnstatestructurefile``   specifies the begin structurefile with bond orders since babel sometimes cannot guess correct bond orders from XYZ file. Default has no value.


``endstatestructurefile``   specifies the end structurefile with bond orders since babel sometimes cannot guess correct bond orders from XYZ file. Default has no value.


``mutlambdascheme``     specifies the lambda scheme (from 0 to 1) of mutation of any parameters between begin and end state keyfiles. Default has no value.


``compareparameters``   specifies to compare parameters with analyze and print parameters if energy difference is > energy cutoff. Default is False.

``energycutoff``    specifies the energy cutoff for comparing parameters. Default is .5.


``rotateframes``    specifies to rotate multipole frames from bgnstatekey to endstatekey. Default is False.


