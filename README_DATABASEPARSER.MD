## Databaseparser README
<ol>
<li>	Determine all single atoms, bonds, angles and torsions from input molecule </li>
<li>	If input forcefield is AMOEBA+ </li>
<ol>
    <li>	Generate map of SMARTS to atom class and atom class to comments from input AMOEBA+ file </li>
    <li>	Match all SMARTS to atom indices in molecule </li>
    <li>	Generate map of atom class to parameters from input file </li>
    <li>	Find best SMARTS match out of all possible for charge flux bond/angle parameters </li>
<ol>
        <li>	For each SMARTS, generate a fragment molecule and use tanimoto coefficient between input molecule and each SMARTS fragment </li>
        <li>	Grab the max tanimoto SMARTS fragment (most similar), if multiple have same value then look at the length of SMARTS string, the longest SMART string (if bond/angle then count length of each SMARTS in tuple) has the most information and this SMARTS is defined as the best match </li>
</ol>
    <li>	Find best SMARTS match out of all possible for bond/angle parameters </li>
    <li>	Construct parameter lines for bond/angle charge flux </li>
    <li>	Read AMOEBA+ non-bonded input files </li>
<ol>
        <li>	Charge penetration </li>
        <li>	Charge transfer </li>
        <li>	vdw </li>
</ol>
    <li>	Match all SMARTS to atom indices in molecule </li>
    <li>	Generate map of atom class to parameters from input file </li>
    <li>	Find best SMARTS match out of all possible for non-bonded parameters </li>
    <li>	Construct parameter lines for non-bonded parameters </li>
</ol>
<li>	Read input SMARTS to solute radii for implicit solvent </li>
<li>	Match all SMARTS to atom indices in molecule </li>
<li>	Find best SMARTS match to solute parameters </li>
<li>	Construct solute parameter lines </li>
<li>	Read database of polarize parameters </li>
<ol>
    <li>	Generate SMARTS to atom class maps </li>
</ol>
<li>	Read database of polarize parameters </li>
<li>	Find best SMARTS match to polarize parameters </li>
<li>	Constructs polarize parameter lines </li>
<li>	Read database of parameters for bond, angle, strbnd, opbend </li>
<ol>
    <li>	Generate SMARTS to atom class maps </li>
</ol>
<li>	Read database of parameters for bond, angle, strbnd, opbend </li>
<li>	Match all SMARTS to atom indices in molecule (including bond, angle, strbnd, opbend) </li>
<li>	Find best SMARTS match for bond, angle, strbnd, opbend </li>
<li>	Construct parameter lines for bond, angle, strbnd, opbend </li>
<li>	Generate list of potential SMARTS for atom, bonds, angles, torsion by iteratively generating all possible fragment indices and consecutively connected atom indices in molecule, then converting to fragment and generating a SMARTS string </li>
<li>	Read file of SMARTS-atom order to tinker element-type description </li>
<li>	Read amoeba09.prm file and map tinker element-type description to tinker type and tinker class </li>
<li>	Match list of possible SMARTS to smarts from database file containing element – tinker type descriptions </li>
<li>	Read external database of SMARTS for bond, angle, strbnd, torsion, opbend, vdw, tortor </li>
<li>	Match all SMARTS in external database to input molecule and return indices that match to SMARTS (for bond, angle, strbnd, torsion, opbend, vdw, tortor) </li>
<li>	For atom (vdw), bond, angle, strbnd, and torsion compare SMARTS matches from external database to matches from iteratively fragmenting input molecule </li>
    <li>	Keep which ever match has the most atoms (longest SMARTS match), unless ~ or * wild cards are detected in SMARTS in which case its possible smaller match is kept </li>
<li>	Generate maps of atom indices to atom type and atom indices to atom class for atom (vdw), bond, angle, torsion (strbnd same as angle) </li>
<li>	 Remove bond and angle matches from tinker element-type description that matched to bond, angle, strbnd, opbend database. </li>
<li>	Find all consecutive rotatable bonds </li>
<li>	Find all missing tortors </li>
<ol>
    <li>	Input tortor external SMARTS matching indices </li>
    <li>	If all 5 atom indices + neighbors exist in SMARTS then tortor is defined as not missing, otherwise it is missing </li>
</ol>
<li>	Find all missing torsion </li>
<ol>
    <li>	This defines which torsions are zeroed out prior to torsion fitting </li>
    <li>	Iterate over torsion indices, SMARTS strings </li>
    <li>	If b-c (a-b-c-d) is in onlyrotbndslist, and if a or d is a hydrogen, if there is an example of a possible a and d that is heavy (non-hydrogen) then skip this torsion (try to transfer hydrogen torsion), else then add to missing torsions list. So if all hydrogen torsion on one side or both sides then fit this torsion. </li>
    <li>	If bond order of b-c is not 1, then skip torsion </li>
    <li>	If the angle between a-b-c and b-c-d is near 180 degrees, then skip torsion </li>
    <li>	If linear, then dihedral angle is not defined, minimize will not work </li>
    <li>	If b-c is in partial double bonds then skip torsion </li>
    <li>	If there are no aromatic atoms (such as if torsion was in a ring) </li>
<li>	If a or d is hydrogen but there is example of a and d as heavy atoms than skip torsion (try to transfer hydrogen) else if no heavy atoms both sides don’t skip torsion. So, if all hydrogen torsion on one side or both sides then fit this torsion. </li>
    <li>	If ~ or * in SMARTS, and no aromatic atoms in torsion (such as if in ring), then this is a poor match and then try to derive parameters for this torsion. So aromatic torsions can transfer ~ * in SMARTS but non-aromatic will decide to derive parameters for. </li>
    <li>	If missing any neighbors of a-b-c-d in SMARTS and not in an aromatic ring, then derive torsion parameters for this torsion. </li>
</ol>
<li>	Find adjacent missing tortors </li>
<ol>
<li>	Use missing torsions as input, if a-b-c-d is missing but b-c-d-e is not missing, then call a-b-c-d-e missing tortor </li>
</ol>
<li>	Find missing vdw parameters </li>
<ol>
    <li>	If SMARTS match does not have atom index + neighbors in SMARTS match then defined as missing vdw parameter </li>
    <li>	Or if ~, * in SMARTS defined as missing parameters </li>
</ol>

<li>	Find examples of atomindices for each missing vdw type </li>
<li>	Find missing bond/angle parameters </li>
<ol>
    <li>	If SMARTS match does not have atom index + neighbors in SMARTS match then defined as missing vdw parameter </li>
    <li>	Or if ~, * in SMARTS defined as missing parameters </li>
</ol>
<li>	Find missing bond/angle parameters </li>
<li>	For all atom, bond, angle, torsion, opbend tinker classes, map tinker classes to poltype classes </li>
<li>	Search amoeba09.prm file for matches to element-type description database </li>
<li>	Correct torsion parameters with pitor energy </li>
<ol>
    <li>	Distribute total pitor energy to all torsion around rotatable bond (the second fold torsion parameter)</li> 
</ol>
<li>	Correct torsion parameters with pitor energy </li>
<li>	Remove high energy hydrogen torsion parameters that dont have any atoms in torsion on a ring and also dont include any hydrogen torsion around partial double bond. </li>
<ol>
    <li>	>= 1 kcal/mol </li>
</ol>
<li>	For angle, bond, opbend, strbnd matches to amoeba09, remove matches that had matches to the bond, angle, strbnd, opbend external SMARTS database. Then add parameters from bond, angle, strbnd, torsion, opbend, vdw, tortor. </li>
<ol>
    <li>	External SMARTS for bond, angle, strbnd, torsion, opbend, vdw, tortor has highest priority (think vdw, torsion, tortor from poltype will be most descriptive SMARTS) </li>
    <li>	Then bond, opbend, strbnd, angle has next highest priority </li>
    <li>	Then amoeba9 matches is last priority </li>
</ol>
<li>	Zero out missing strbnd parameters </li>
<li>	Assign bond/angle guesses for missing angle bond parameters </li>
<ol>
    <li>	Using element and valence to determine guesses </li>
    <li>	Taken from valence.f tinker program </li>
</ol>
<li>	Modify angle keywords that need anglep (planar) keyword </li>
<li>	Modify angle and bond parameters to use QM optimized bond and angle lengths </li>
<li>	Add missing hydrogen torsions that were screened due to high energy transfers to list of torsion that needs to be zeroed out </li>
<li>	Zero out all torsion parameters that are defined as missing </li>
<li>	For missing aromatic torsions, default to benzene parameters </li>
<li>	Construct comment lines for all parameter lines </li>
<li>	Write out lists of missing torsions, missing vdw, missing tortors </li>
</ol>
