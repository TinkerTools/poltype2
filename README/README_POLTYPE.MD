
## Poltype Protocol

* See original publication https://doi.org/10.1007/s00214-012-1138-6
* Below is a high level flow diagram for Poltype
* Vdw fitting is not yet implemented
* For more details about individual protocols in flow diagram, see the corresponding README files. 
![Poltype Protocol](Images/PoltypeProtocol.PNG)
* A more detailed representation of filenames is presented below
![Poltype Protocol Detailed](Images/PoltypeProtocolDetailed.PNG)

<ol>
<li>	Initialize defaults for keywords in the class PolarizableTyper </li>
<li>	Read input keywords from poltype.ini file </li>
<li>	Change some QM method keywords depending on if using gaussian/psi4 (sometimes DFT functional keywords are not the same between two programs) </li>
<li>	If tinker executables are not in PATH variable, then add .x extension to each of them (sometimes users compile them differently) </li>
<li>	Print copyright information </li>
<li>	Initialize paths need for all executables, check tinker version </li>
<ol>
    <li>Search for g03, g09 and g16 </li>
    <li>Search for cubegen executable </li>
    <li>Search for formchk executable </li>
    <li>If tinker version is not > default (8.8.3) then crash with error </li>
    <li>Search for TINKERDIR variable (else program executable should be in PATH variable) for poledit, potential, minimize, analyze and superpose programs </li>
    <li>Search for GDMA executable and crash if using gaussian for SP computations and cannot find executable </li>
<ol>
        <li>	If using psi4 then GDMA is called from psi4 program instead of executable </li>
</ol>
    <li>Search for psi4 in PATH and if cannot find then crash </li>
    <li>Search for GAUSS_SCRDIR in environment and if cannot find while using gaussian for SP or OPT then crash </li>
</ol>
<li>	Initialize filenames to be used </li>
<li>	If amoebabio18.prm path and modified protein pdb and unmodified protein pdb are given as inputs then	</li>
<ol>
    <li>Generate poltype input molecule (including parts of backbone) </li>
    <li>Generate parameters </li>
    <li>Generate modified protein xyz and key </li>
</ol>
<li>	Generate parameters </li>
<ol>
    <li>Check if input processors is more than default of 8 (too many wont give any extra benefit) </li>
    <li>Initialize openbabel and rdkit mol object and read in input structure file </li>
    <li>Create scratch directories for this molecule if they don’t already exist </li>
    <li>Generate the extended conformer (some molecules structures are floppy, and start compact, want to find most extended for multipoles to not pick up any extra intermolecular interactions). </li>
<ol>
        <li>	Use rdkit to generate many confomers (via some stochastic/genetic algorithm) </li>
        <li>	For each conformer then compute longest distance between two atoms </li>
        <li>	Pick conformer with maximum longest distance </li>
        <li>	Rdkit sometimes gives incorrect bond lengths (via some cheap minimization scheme), so just readjust bond lengths back to the original bond lengths via vector algebra </li>
</ol>
    <li>	If input is given for keyfile to be added to parameter database, then add parameter lines to amoeba09.prm file and terminate poltype </li>
    <li>	If Iodine is detected in the molecule and Gaussian is detected in the environment, then use gaussian for SP and OPT computations. </li>
    <li>	If total charge is 0 but there is an atom with formal charge, then assume this is zwitterion (need better detection). </li>
<ol>
        <li>	Then if Gaussian is detected, use Gaussian for OPT and SP (psi4 has poor PCM optimizer compared to Gaussian). Performance may only be significant for larger molecules. </li>
        <li>	Change QM method names again (must change when changing from psi4 to gaussian and vice versa) </li>
</ol>
    <li>If the number of atoms is less than 25, then turn the fragmenter module off </li>
    <li>Delete any cartesian XYZ files (if user leaves these in input folder, will make tinker confused) </li>
    <li>Determine the symmetry classes via symmetry.py </li>
    <li>Try to optimize geometry (in optimization.py), if fail then try again with more opt cycles </li>
<ol>
        <li>	Start with 5 cycles </li>
        <li>	After each completion of 5 cycles check the bond lengths compared to input structure. If any bond deviates more than .2 angstroms, then optimization begins again with +2 more optimization cycles until deviations are less than the tolerance </li>
</ol>
    <li>	Find partial double bonds (defined as bonds in amide/acids) via SMART string search, then exclude these from torsion parameter derivation (transfer and don’t fit). </li>
    <li>	If .key_4 is missing, then call databaseparser.py and using optimized structure as input, grab bond, angle, strbnd, torsion and vdw parameters from database files. </li>
<ol>
        <li>	Determine missing torsions </li>
        <li>	Determine missing tor-tors </li>
        <li>	Determine missing vdw parameters </li>
        <li>	Grab torsion parameter guesses (will be used for non-aromatic ring puckering if keyword is turned on). </li>
</ol>
    <li>Compute DMA SP in electrostaticpotential.py </li>
    <li>Compute GDMA (or read from psi4 output file) and create .gdmaout file in electrostaticpotential.py </li>
    <li>Generate poledit inputfile (.peditin) in multipole.py </li>
    <li>Call poledit program (inputs .gdmaout and .peditin) </li>
    <li>Generates .xyz, .key and .grid file </li>
    <li>Add header to key file </li>
<ol>
        <li>	Turn all parameters off except multipole and polarize </li>
        <li>	Add parameter line (such as amoebabio18.prm) </li>
</ol>
    <li>Compute ESP SP in electrostaticpotential.py </li>
    <li>Generate ESP .pot file for tinker (in electrostaticpotential.py) </li>
    <li>Used for fitting </li>
    <li>Average multipoles via same symmetry class (in multipole.py) </li>
<ol>
        <li>	Generates .xyz_2 and .key_2 </li>
</ol>
    <li>If .key_3 does not exist then call electrostatic potential fitting (in electrostaticpotential.py) </li>
<ol>
        <li>	Generates .key_3 </li>
</ol>
    <li>Dump parameters found from databaseparser.py into key_4 (using key_3 as input) (in databaseparser.py) </li>
    <li>If PCM was turned on, then add solvate GK keyword to key_4 file </li>
    <li>Using missing torsions as input, determine which rotatable bonds will have torsion parameters derived for (in torsiongenerator.py) </li>
    <li>Remove duplicate rotatable bonds (via symmetry classes) (in torsiongenerator.py) </li>
    <li>Initialize all variable torsions (those to not be restrained during minimization) as no torsions (empty list). </li>
    <li>Prepare rotatable bonds in input format for fragmenter </li>
    <li>If vdw keyword is turned on then take missing vdw atom indices and prepare them in input format for fragmenter </li>
    <li>Determine the angle increment needed for each torsion (12 degrees etc..) </li>
    <li>If tortor keyword is on find all consecutive rotatable torsions (for bond a-b there are many torsions and there are consecutive torsions around b-c) (in torsiongenerator.py) </li>
<ol>
        <li>	Prepare inputs of tor-tors for fragmenter (in torsiongenerator.py) </li>
</ol>
    <li> If refinearoringtors keyword is on, then iterate over all non-aromatic rings and provide inputs for fragmenter module (in rings.py) </li>
    <li>Generate WBO index matrix for parent (in fragmenter.py) </li>
    <li>Draw a 2D image of molecule and bonds that are missing parameters (in fragmenter.py) </li>
    <li>Generate fragments via iteratively comparing WBO parent indices to fragment indices ( in fragmenter.py) to reach some threshold difference in WBO value between fragment and parent. </li>
    <li>Spawn child poltype jobs for each fragment (in fragmenter.py) </li>
    <li>Grab vdw, torsion and tortor parameters from fragments </li>
    <li>If a fragment or if fragmenter turned off then Generate torsion parameters ( in torsiongenerator.py) </li>
    <li>If a fragment or if fragmenter turned off then fit torsion parameters from QM-MM energy profiles (in torsionfit.py) </li>
    <li>If vdwfitting keyword is turned on and is a fragment then derive vdw parameters from QM-MM energy profiles for missing vdw atom indices (in vdwfit.py) </li>
<ol>
        <li>	Also if fragmenting is turned off but vdwfitting turned on then derive vdw parameters for missing vdw atom indices </li>
</ol>
    <li>Check torsion parameters transferred from fragments to see if any have all 0’s </li>
<ol>
        <li>	If so then crash program (this was used early on when fragmenter was buggy) </li>
</ol>
    <li>Minimize final MM structure then compare with QM (in optimization.py) </li>
<ol>
        <li>	Check to see if RMSD of structure is <1 else crash program </li>
</ol>
    <li>If PCM was turned on remove solvate GK from key file </li>
    <li>Compare MM and QM dipole moments (in electrostaticpotential.py) </li>
    <li>Delete scratch files </li>
    <li>Write out termination signal to poltype log </li>
</ol>
</ol>
