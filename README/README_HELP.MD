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

``addhydrogentocharged`` specifies to add hydrogens to carbon with valence less than 4 and nitrogen with valence 2. Default is True.

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

``optpcm`` specifies to use PCM for geometry optimization. [True|False|auto] If set to auto, PCM will be used for zwitterion. Default is auto.

``pcm_auto`` specifies to automatically activate PCM for zwitterions. If False, PCM will not be automatically activated. Default is True.

``dontusepcm`` equivalent to `pcm_auto=False`.

``use_gausoptonly``  specifies to only use gaussian for initial geometry optimization and Psi4 for all other QM steps. Default is False. 

``generateextendedconf`` specifies to use extended conformation for geometry optimization with torsions restrained. Default is True.

``userconformation`` specifies to use the user input conformation for geometry optimization with torsions restrained. Default is False.

``userxyzgeometry`` specifies an xyz file to use user input geometry (usually already optimized), Poltype will skip QM opt and continue with DMA step. Default is ''.

``sp2aniline`` specifies to generate flat geometry of c-NH2 (only works with `use_gaus` keyword). If True, poltype will also assign `anglep` and `opbend` parameters. Default is True.

``optloose`` specifies to use GAU_LOOSE optimization convergence criteria for psi4 optimization (both first geometry optimzation and for torsion geometry optimization)

``optconvergence`` specifies optimization convergence. LOOSE, GAU_LOOSE criteria (identical to the keyword ``optloose``); VERYLOOSE, only using gradient thresholds from GAU_LOOSE criteria. Default is LOOSE. 

#### Distributed Multipole Analysis and Electrostatic Single Point

``new_gdma`` specifies to use the new grid-based DMA algorithm (experimental feature). Default is False.

``gdmacommand_*`` add/modify GDMA commands. For example, `gdmacommand_Radius_S=0.80` will add the command `Radius S 0.80` to GDMA input file.

``usepoleditframes`` specifies to use POLEDIT frames rather than poltype frame detection. Default is True.


``dmamethod``  specifies which quantum method to use for distributed multipole analysis (DMA). Default is MP2.
                  
``dmabasisset`` specifies basis set to use for quantum that generates initial multipoles. Default is 6-311G**.         
                  
``espmethod``   specifies which quantum method to use for electrostatic fitting. Default is MP2.

``espbasisset`` specifies basis set to use for electrostatic fitting. Default is aug-cc-pvtz.

``sameleveldmaesp`` specifies to use only do QM once (DMA) and copy files for ESP. Default is False.

``adaptiveespbasisset`` specifies to downgrade ESP basis set from aug-cc-pvtz to aug-cc-pvdz when molecule have >=20 atoms. Default is False.

``numespconfs`` specifies number of conformations to be used for ESP potential fitting. Rdkit conformors are evaluated via MMFF94 energies and numespconfs-1 most favorable conformations are picked to add to fitting (with first being most extended conformation). Default is 1.

``espextraconflist`` specifies list of conformors to add to extended conformation for ESP fitting. Default is to just use extended conformation for ESP fitting. Files can be SDF or any format as long as atom order is same as input molecule atom order.

``indextompoleframefile`` specifies file of multipole frame definitions for each line (exactly how frames are defined in first key file). This will be used instead of poltype frame detection if specified.

``scalebigmultipole``   when big multipoles seen in DMA, specifies to scale down, then to fit high level ESP. Default is False.

``fragbigmultipole``   when big multipoles seen in DMA, specifies to use a smaller fragment to derive multipole, then to fit high level ESP. Default is True.

``atomidsfordmafrag``   specifies a list of atoms to be thought as "BID multipole atoms" to use `fragbigmultipole`. Default is [].


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

``prmmodfile`` specifies (comma-separated list of) parameter modification patch file(s). See `ParameterFiles/dma4_hfe2023.mod` for example.

#### Torsion Parameterization

``xtbtorresconstant`` specifies torsion restraint force constant for restraining dihedrals in xtb. Default is 5 kcal/mol/deg^2.

``rotalltors``  specifies poltype to spin all rotatable bonds instead of transferring torsion parameters if found in database. Default is False.
                  
``torsppcm``    specifies poltype to use PCM for torsion SP. [True|False|auto] If set to auto, PCM will be used for zwitterion. Default is auto. 

``toroptpcm``   specifies poltype to use PCM for torsion optimization. [True|False|auto] If set to auto, PCM will be used for zwitterion. Default is auto. 
                  
``dontdotor``   specifies poltype to skip torsion parameterization. Default is False.
                  

``torsionrestraint``  specifies how many kcals/mol/rad^2 to restrain other torsions around rotatable bond to be at current dihedral angle. Default is 1641.4 (.5 kcal/mol/deg^2).

``foldnum``           specifies the number of fourier terms used in torsion fitting. Default is 6.

``tordatapointsnum``  specifies the number of quantum jobs to sample along the dihedral angle for all torsions. Default is 12 if number of torsion parameters to fit is less than 12, else automatically adds more points.

``toroptmethod`` specifies the opt method used for torsion parameterization. Default is xtb. ``toroptmethod=xtb`` will use semi empirical qm methods (not psi4 or gaussian). ``toroptmethod=ANI`` will use ANI for structure optimization. If ``toroptmethod=AMOEBA`` then AMOEBA minimized structure will be used instead of QM opt structure.

``toroptmethodlist`` specifies a list of opt methods to use and can compare QM energy torsion profiles. By default ``toroptmethodlist`` only has ``toroptmethod`` as just one item in list.

``torspmethodlist`` specifies a list of opt methods to use and can compare QM energy torsion profiles. By default ``torspmethodlist`` only has ``torspmethod`` as just one item in list.

``xtbmethod`` specifies which xtb semi-empirical method to use, must have ``xtb`` turned on for either toroptmethod or torspmethod. Default is 2 (this corresponds to GFN2-xTB). 

``toroptbasisset``    specifies which basis set to use for geometry optimization during torsion parameterization. Default is 6-31G*. 

``torspmethod`` specifies the SP method for torsion parameterization. Default is wB97XD.``torspmethod=xtb`` will use semi empirical qm methods (not psi4 or gaussian). ``torspmethod=ANI`` will use ANI for SP energies.

``torspbasisset``  specifies which basis set to use for SP during torsion parameterization. Default is 6-311+G*.
                  

``onlyrotbndslist``   specifies to only rotate bonds in the list, onlyrotbndlist=1 2, 3 4, 5 6 . Default has no value.


``dontrotbndslist``   specifies to not rotate bonds in the list, onlyrotbndlist=1 2, 3 4, 5 6 . Useful for modified AA lib, to tell poltype not to rotate backbone dihedrals. Default has no value.


``tortor``            specifies to sample the torsion-torsion potential energy surface for two adjecent and rotatable dihedrals. Default is False.

``onlyrottortorlist`` specifies list of tortors to rotate. "onlyrottortorlist=b1 c1 d1 , b2 c2 d2" where b1,c1,d1 belong tor tortor a1,b1,c1,d1,e1. Default has no value.
  

``refinenonaroringtors`` specifies to pucker and refine torsion parameters on non-aromatic rings. Default is False.


``nonaroringtor1Dscan`` specifies to cut a non-aromatic ring and spin the torsion for a 1D scan. Default behaviour is to transfer alkane torsion parameters.
                    
``fitfirsttorsionfoldphase`` specifies to add a phase fitting parameter to the first cosine term in the torsion fitting                           function. Default is False.

``maxtorresnitrogen`` specifies number of dihedral restraints per rotatable bond for rotatable bonds that contain nitrogen dihedrals with hydrogens. Default is 1. For comparing different optimizers that may shift heavy atoms differently can set higher.  

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

