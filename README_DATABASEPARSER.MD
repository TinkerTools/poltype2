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
    <li>	If bond order of b-c is not 1 and bond is not in a ring, then skip torsion </li>
    <li>	If the angle between a-b-c and b-c-d is near 180 degrees, then skip torsion </li>
    <li>	If linear, then dihedral angle is not defined, minimize will not work </li>
    <li>	If b-c is in partial double bonds then skip torsion </li>
    <li>	If missing any neighbors of a-b-c-d in SMARTS string or *, ~ in SMARTS string, then determine this to be poorly transferred parameters. </li>
    <li>	If ring bond and missing torsion from database: </li>
<ol>
             <li> If a-b-c-d all non-aromatic  </li>
<ol>
                  <li> If a or d is H and and there exist other a/d such that a,d is not H, then just transfer whatever matches for H torsion </li>
<ol>
                          <li> Exception when transferring high energy H torsion, then zero out and fit  </li>
</ol>
                  <li> Else, chop ring and derive parameters  </li>
</ol>
             <li> b-c aromatic (or b-c have hybridization of 2). Below, lowercase leter means aromatic atom or ring atom with hyb of 2. Letter x refers to non-hydrogen atom.  </li>
<ol>
                    <li> If H-x-x-x transfer from H-c-c-c  </li>
                    <li> If H-x-x-H transfer from H-c-c-H  </li>
                    <li> If x-x-x-x transfer from c-c-c-c  </li>
</ol>
             <li> If b aromatic and c not aromatic or (b non-aromatic and c aromatic).   </li>
<ol>
                    <li> If H-X-X-X transfer from H-C-C-C  </li>
                    <li> If H-X-X-H transfer from H-C-C-H  </li>
                    <li> If X-X-X-X transfer from C-C-C-C  </li>
</ol>
             <li> If b and c non-aromatic but a or d is aromatic   </li>
<ol>
                    <li> If H-X-X-X transfer from H-C-C-C  </li>
                    <li> If H-X-X-H transfer from H-C-C-H  </li>
                    <li> If X-X-X-X transfer from C-C-C-C  </li>
</ol>
</ol>
             <li> If NOT ring bond and missing torsion from database:   </li>
<ol>
                  <li> If a or d is H and and there exist other a/d such that a,d is not H, then just transfer whatever matches for H torsion </li>
<ol>
                          <li> Exception when transferring high energy H torsion, then zero out and fit  </li>
</ol>
                  <li> Else, derive parameters for torsion  </li>
</ol>
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

### Example keyfiles comments
*  ['[#16]-[H]', [1]] this first SMARTS, describes a fragment of the input molecule containing the vdw atom type and a neighbor. This SMARTS is the closest match found to all amoeba0 SMARTS matches. The 1 indicates first atom in SMARTS string. Similarily, the [2] in the second SMARTS string (from amoeba09 molecule), indicates the second atom in that SMART string is where the parameter is transferred from. [('S', '"Ethyl Sulfide S"')] this shows the tinker type description from amoeba09 atom.
* "Missing vdw parameters, will attempt to fit parameters", this means that the SMARTS '[#16]-[H]' does not contain all neighbors of the sulfur atom in the input molecule and thus the parameters are determined to be missing and not transferrible enough, so parameters will be derived. 

```
# matching SMARTS from molecule  ['[#16]-[H]', [1]] to SMARTS from parameter file ['[H]-[#16]-[#6](-[#6](-[H])(-[H])-[H])(-[H])-[H]', [2]] with tinker type descriptions [('S', '"Ethyl Sulfide S"')]
# Missing vdw parameters, will attempt to fit parameters
vdw 401 4.0050 0.3550

```
* "Missing bond parameters, assigning default parameters", this means that bond parameters SMARTS "[#16]-[#6]", did not contain all neighbors in input molecule and so the parameters are determined not to be transferrible and default parameters based on element and valence are assigned.

```
# matching SMARTS from molecule  ['[#16]-[#6]', [1, 2]] to SMARTS from parameter file ['[#16](-[#8-])(=[#8])(=[#8])-[#6](-[#6](-[#6](-[H])(-[H])-[H])(-[H])-[H])(-[H])-[H]', [1, 5]] with tinker type descriptions [('S', '"Ethyl Sulfonate SO3-"'), ('C', '"Ethyl Sulfonate CH2"')]
# Missing bond parameters, assigning default parameters
bond         401   402          250.0     1.8015
```
* "WARNING WILDCARDS USED IN SMARTS PARAMETER MATCHING", another way that parameters are determined to be non-transferable and default parameters involving element and valence are assigned.
```
# matching SMARTS from molecule  ['[#6]~[#6]~[#16]~[H]', [1, 2, 3]] to SMARTS from parameter file ['[H]-[#16]-[#6](-[#6](-[H])(-[H])-[H])(-[H])-[H]', [4, 3, 2]] with tinker type descriptions [('C', '"Ethyl Sulfide CH3"'), ('C', '"Ethyl Sulfide CH2"'), ('S', '"Ethyl Sulfide S"')]
# WARNING WILDCARDS USED IN SMARTS PARAMETER MATCHING
# Missing angle parameters, assigning default parameters
angle       403   402   401      53.2     119.02
```
* New database for bond,angle,strbnd and opbend supersedes the amoeba09 database matches and has comments such as " updated valence parameter database match, comments=car, car, car, SMARTS match = [c] [c] [c]"
```
# updated valence parameter database match, comments=car, car, car, SMARTS match = [c] [c] [c]
angle 404 403 402 74.630713 118.66
```
* When parameters are missing in the database for strbnd, they are removed and you will see a comment such as "Missing strbnd parameters, zeroing out parameters"

```
# matching SMARTS from molecule  ['[H]-[#16]-[#6]', [1, 2, 3]] to SMARTS from parameter file ['[H]-[#16]-[#6](-[#6](-[H])(-[H])-[H])(-[H])-[H]', [1, 2, 3]] with tinker type descriptions [('H', '"Ethyl Sulfide HS"'), ('S', '"Ethyl Sulfide S"'), ('C', '"Ethyl Sulfide CH2"')]
# Missing strbnd parameters, zeroing out parameters
strbnd       407   401   402       0       0
```
* If missing opbend parameters from database, MM3 forcefield default parameters are used and a comment in the key file will say "WARNING DEFAULT MM3 OPBEND VALUES USED "
# matching SMARTS from molecule  ['[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1', [1, 2]] to SMARTS from parameter file ['[#8](-[#6]1:[#6](:[#6](:[#6](:[#6](:[#6]:1-[H])-[H])-[H])-[H])-[H])-[H]', [2, 3]] with tinker type descriptions [('C', '"Phenol C1-OH"'), ('C', '"Phenol C2"')]
# WARNING DEFAULT MM3 OPBEND VALUES USED 
opbend 402 403 0 0 14.39
```
* If missing aromatic torsion (or b-c bond in a-b-c-d both are in ring with hybridization of 2) from database benzene parameters are transferred 
```
# Transferring from ['Benzene HC', 'Benzene C', 'Benzene C', 'Benzene C']
torsion 409 404 405 404 0.55 0.0 1 6.187 180.0 2 -0.55 0.0 3
```
* If missing torsion parameters from database, a comment will appear saying "Missing torsion parameters, will attempt to fit parameters"
 
```
# matching SMARTS from molecule  ['[#8]~[#16]~[#6]~[#6]~[#6]~[H]', [1, 2, 3, 4]] to SMARTS from parameter file ['[#16](-[#8-])(=[#8])(=[#8])-[#6](-[#6](-[#6](-[H])(-[H])-[H])(-[H])-[H])(-[H])-[H]', [2, 1, 5, 6]] with tinker type descriptions [('O', '"Ethyl Sulfonate SO3-"'), ('S', '"Ethyl Sulfonate SO3-"'), ('C', '"Ethyl Sulfonate CH2"'), ('C', '"Propyl Sulfonate Me-CH2"')]
# WARNING WILDCARDS USED IN SMARTS PARAMETER MATCHING
Missing torsion parameters, will attempt to fit parameters
torsion 406 401 402 403 0 0.0 1 0 180.0 2 0 0.0 3

```
