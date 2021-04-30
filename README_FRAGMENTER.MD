## Fragmenter README
* Quantum computations scale exponentially with increasing basis set (number of electrons) in your system
* Some large molecules will have very poor torsion fits because there are lots of forces from other parts 
  of the molecule (steric clashes). 
* It is better to derive the "intrinsic" torsion of rotatable bonds without steric clashes is to fragment the   molecule. However we want to ensure that the electron density in the parent molecule nearby the inner most    rotatable bond (b-c in torsion a-b-c-d) is similar to the electron density of the fragment for the same bond

* The fragmenter library implements the Wiberg bond order to have estimates of electron density on rotatable    bonds of interest. The relative Wiberg bond order (WBO) between parent and fragment is used as a cost         function to obtain the best fragment that is smaller than the parent (has relative WBO less than some tolerance)

* A flow diagram is presented below

* parenthesis for atom (e) is if user chooses to fragment with two adjacent torsions rather than just one, where for the atoms a-b-c-d-e, the rotatable bonds would be b-c and c-d. Currently fragmenting with two adjcent rosions is not yet supported. 

![Fragmentation Flow Diagram](Images/FragmentationProtocol.PNG)


* An example below identifies missing torsions by highlighting corresponding rotatable bonds. An example of the fragmenting scheme is shown for one of the missing torsions below. A tolerance of .001 was used.

* If two missing torsions from the parent have the same converged fragment, then only one is chosen to submit the poltype job and both rotatable bonds are spun respectively for that poltype job.

* Bond topology need to be correct for formal charge assignment on fragments in order to get the total charge of the fragment correctly. 

![Example Molecule Image](Images/ml188.png)

![Fragment Growth Example Image](Images/FragmentGrowthWithIndex_Bnd_11-8_Index_0.png)

* FragmentGrowthWithIndex_Bnd_b_c_Index_0.png  Tracks the growth of fragment along while attempting to                                                       minimize differences in Wiberg Bond Order between each fragment                                               and the parent until the tolerance threshold is reached. If                                                   WithIndex in name, this means it will show atoms with their                                                   index numbers. If WithoutIndex in name, it will show atoms                                                    without their index numbers. The last molecule in image is                                                    the molecule that will be used to derive torsion paramaters                                                   for the parent molecule.



![Combinations Example Image](Images/CombinationsWithIndex_Bnd_11-8_Index_0.png)

* CombintationsWithIndex_Bnd_b_c_Index_0.png   There are different combinations of ways to grow a fragment out                                               and each way is attempted before choosing the minimum                                                         difference in Wiberg bond order between parent bond and                                                       fragment bond. This minimum fragment is added to the                                                          FragmentGrowth image (last fragment on Combinations image).


* Several dictionaries are available in the fragmenter directory for troubleshooting
1) parentindextofragindex.txt        map from parent atom indices to fragment atom indices
2) classkeytosmartsposarray.txt      map from parents torsion key (four type numbers) to the four indices of the torsion atoms in the SMARTS string.
3) fragclasskeytoparentclasskey.txt  map from fragment torsion class key to parent torsion class key
4) parentsymclasstofragsymclass.txt  map from parent type class to fragment type class
5) classkeytotorsionindexes.txt      map from parent torsion key to the torsion atom indices
### Protocol
<ol>
<li>	Generate Fragments </li>
<ol>
    <li>	For each input list of atom indices (could be for torsion, tortor or vdw or non-aromatic ring atom indices) </li>
    <li>	Iterate over input list of atom indices </li>
    <li>	If not aromatic ring atom indices then </li>
<ol>
        <li>	Generate fist pass atom indices via adding neighbors of atom indices in parent from input indices </li>
</ol>
    <li>	If vdw fragment or non-aromatic ring indices fragment, then set fragment growth cycles to 0 (don’t grow out) </li>
    <li>	Generate initial fragment from atom indices </li>
    <li>	Generate WBO matrix for fragment via HF/MINIX </li>
    <li>	Compute the difference between parent WBO and fragment WBO around rotatable bond and save that </li>
    <li>	Draw 2D images and highlight bonds of interest with absolute and relative WBO values </li>
    <li>	If maxgrowthcycles is not 0, then try growing the fragment out </li>
    <li>	Find equivalent fragments from all fragments generated </li>
    <li>	Copy first equivalent fragment structure to other equivalent structures (sometimes babel has same molecule but different order of atom indices) </li>
</ol>
<li>	Grow fragment out </li>
<ol>
    <li>	Generate all possible sets of newly grown atom indices </li>
<ol>
        <li>	If multiple bonds cut then many directions to grow out and thus many combinations of ways to grow out </li>
</ol>
    <li>	For each possible fragment indices , generate the fragment and compute the WBO matrix, compute difference between parent and fragment WBO values and save it. </li>
    <li>	Draw 2D images </li>
    <li>	Find the fragment with the smallest difference between parent and fragment WBO values and use this fragment as the final fragment for this iteration. </li>
    <li>	Keep iterating until either maxgrowth cycles (4) are reached or threshold relative WBO is less than or equal to the tolerance (.05). </li>
</ol>
<li>	Generate Fragment </li>
<ol>
    <li>	Take input of parent atom indices. </li>
    <li>	For each atom index in input, grab parent atom object and add properties to new mol object. </li>
    <li>	Iterate over bonds in parent and if both atom indices are in new fragment, then add that bond to new mol object. </li>
    <li>	Add hydrogens (if made cuts, now add hydrogens to bonds). </li>
<ol>
        <li>	Don’t keep hydrogens added away from where bond was cut (such as if there was a -1 formal charge). </li>
</ol>
    <li>	Generate final parent index to fragment index map. </li>
    <li>	Assign total chare to fragment based on element and bond connectivity. </li>
</ol>
<li>	Find equivalent fragments via checking if fragment SMART strings match each other (and same number of atoms). </li>
<li>	For each equivalent rotatable bond fragment set. </li>
<ol>
    <li>	Use the first fragment in equivalent set of fragments for deriving parameters. </li>
    <li>	Convert parent rotatable bond indices to fragment indices via map from parent index to fragment index. </li>
    <li>	Collect fragment indices and add to poltype.ini file to only tell fragment poltype job to derive parameters for specific rotatable bonds. </li>
    <li>	Generate parent symmetry class to fragment symmetry class map via parent index to fragment index map and map from parent index to symmetry class and fragment index to symmetry class map. </li>
    <li>	If equivalent fragment. </li>
<ol>
        <li>	Grab parent torsions around rotatable bond. </li>
        <li>	Convert parent indices to fragment indices </li>
        <li>	Generate parent class key to fragment class key and save in .txt file. </li>
        <li>	Determine SMARTS, smartsposarray and torsion indexes save in dictionaries and output to .txt files. </li>
</ol>
</ol>
<li>	Set up fragment poltype jobs. </li>
<ol>
    <li>	Take important keywords from parent input file and add to fragment poltype.ini file. </li>
    <li>	Check if parent and fragment are the same molecule, if so copy all QM data from parent to fragment folder and rename files appropriately. </li>
</ol>
<li>	Grab vdw, torsion and tortor parameters from fragments. </li>
<ol>
    <li>	For every set of equivalent rotatable bonds (each fragment is assigned a rotatable bond, which it starts from, some fragments are equivalent but start growing from different rotatable bonds).</li>
<ol>

<li>	For each rotatable bond, if not a vdw fragment, then grab all torsions around rotatable bond from parent molecule and generate torsion and tortor class keys (defined by type numbers for each atom index). Save class keys in an array to search for equivalent fragment class keys. </li>
</ol>
    <li>	For each fragment, navigate to poltype job folder and open dictionaries (stored as .txt files) such as. Index means atom index, class key is a string of symmetry type numbers, smarts in SMARTS string, smartsposarray indicates which atom indices in smart string correspond to the torsion/tortor/vdw atom indices. </li>
<ol>
        <li>	Parentindextofragindex, Parentsymclasstofragsymclasses, Classkeytosmartsposarray, Classkeytosmarts, Parentclasskeytofragclasskey, Classkeytotorsionindexes, Parenttortorclasskeytofragtortorclasskey, Tortorclasskeytosmartsposarray, Tortorclasskeytosmarts, Tortorclasskeytotorsionindexes </li>
        <li>	For every torsion in parent torsions, convert parent class key to fragment class key via dictionaries </li>
        <li>	For every tortor in parent torsions, convert parent class key to fragment class key via dictionaries </li>
        <li>	Open .key_5 file and </li>
<ol>
<li>	If torsion/tortor/vdw is detected in line then check if fragment class key is in torsion line, if it is then </li>
<ol>
    <li>	Determine the SMART string, torsion indexes, smartsposarray  (used for comments in final key file) via above dictionaries and then save parameters and this information for later adding back to parent key file </li>
</ol>
</ol>
</ol>
</ol>
    <li>	Open parent .key_5 file </li>
<ol>
        <li>	If torsion/vdw/tortor in line </li>
<ol>
<li>	Then construct comment line and grab parameters via map from parent to fragment parameters, then construct new parameter line. Write lines to keyfile. </li>
</ol>
</ol>
</ol>
