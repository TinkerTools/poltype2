[Parameterization Output Files](#parameterization-output-files)


* [Final XYZ File](#final-xyz-file)

*  [Final Key File](#final-key-file)

    * [Atom Type Definitions Example](#atom-type-definitions-example)

    * [Van der Waals Parameter Definitions Example](#van-der-waals-parameter-definitions-example)

    * [Bond Parameter Definitions Example](#bond-parameter-definitions-example)

    * [Angle Parameter Definitions Example](#angle-parameter-definitions-example)

    * [Stretch Bend Parameter Definitions Example](#stretch-bend-parameter-definitions-example)

    * [Out of Plane Bend Parameter Definitions Example](#out-of-plane-bend-parameter-definitions-example)

    * [Torsion Parameter Definitions Example](#torsion-parameter-definitions-example)

    * [Solute Parameter Definitions Example](#solute-parameter-definitions-example)

    * [Polarize Parameter Definitions Example](#polarize-parameter-definitions-example)

    * [Multipole Parameter Definitions Example](#multipole-parameter-definitions-example)

* [Poltype Log File](#poltype-log-file)

* [OPENME Plots](#openme-plots)

### Parameterization Output Files


#### Final XYZ File
```
    10
     1  O      1.251722   -1.128767    0.434331      404    5    10
     2  O      0.844886    1.092397    0.237216      406    5
     3  N     -1.815560    0.588102   -0.324599      403    4     8     9
     4  C     -0.849874   -0.492647   -0.449814      401    3     5     6     7
     5  C      0.470594   -0.061298    0.132685      402    1     2     4
     6  H     -1.195550   -1.387129    0.076104      405    4
     7  H     -0.620445   -0.803981   -1.483080      405    4
     8  H     -2.537209    0.496409   -1.034572      407    3
     9  H     -1.333966    1.469038   -0.494689      407    3
    10  H      2.154510   -0.772681    0.564391      408    1
```


* Total atom number is on the first line.
* The first column is the atom index.
* The second column is the atomic symbol.
* The 3 -5th columns are x,y,z coordinates in Angstrom.
* The 6th is the “atom type” defined in the *.key file. This is the index tinker uses to assign parameters from the key/parameter file.
* The 7th – last columns are lists of atom indices that are connected to the current atom index.

#### Final Key File

##### Atom Type Definitions Example
```
atom          404    404    O     "glycine             "         8    15.999    2
atom          406    406    O     "glycine             "         8    15.999    1
atom          403    403    N     "glycine             "         7    14.007    3
atom          401    401    C     "glycine             "         6    12.011    4
atom          402    402    C     "glycine             "         6    12.011    3
atom          405    405    H     "glycine             "         1     1.008    1
atom          407    407    H     "glycine             "         1     1.008    1
atom          408    408    H     "glycine             "         1     1.008    1
```
* First number is the "type" number and the second number is the "class" number. 
* Multipole and Polarize parameters always use type numbers due to the highly specific electrostatic envioronment. 
* All other parameter types use "class" numbers and are less specific.
* By default, poltype uses the same class numbers as type numbers.

##### Van der Waals Parameter Definitions Example

```
# matching SMARTS from molecule  [['[#7](-[#6](-[#6])(-[#1])-[#1])(-[#1])-[#1]', [2]]] to SMARTS from parameter file [['[#7](-[#6](-[#6](-[H])(-[H])-[H])(-[H])-[H])(-[H])-[H]', [2]]] with tinker type descriptions [[('C', '"Ethyl Amine CH2"')]]
# [401] = [[4]]
vdw 401 3.8200 0.1010
```
* All type lines have a line above indicating which indices it corresponds to ([401] = [[4]]), where type number 401 has indices of 4 that correspond to it.
* This type of comment is a match to the amoeba09 database of parameters.
* ['[#7](-[#6](-[#6])(-[#1])-[#1])(-[#1])-[#1]', [2]] the first item in this list is a SMARTS string matching to the input molecule, the second item in the list specifies which atom in order (start counting from 1 on the left) that the match for the vdW atom corresponds to. 
* [['[#7](-[#6](-[#6](-[H])(-[H])-[H])(-[H])-[H])(-[H])-[H]', [2]]] similarly, the first item in this list is a SMILES from a molecule in the amoeba09 database. The seocnd item in the list is the atom in the SMARTS that the match corresponds to.
* [[('C', '"Ethyl Amine CH2"')]] this is a list of the atom class descriptions that are matched from the amoeba09 database
* The first number in the vdW parameter line is radius and the second is the depth parameter


##### Bond Parameter Definitions Example
```
# updated valence parameter database match, comments=C=O, sp2 carbon, carboxylic ester OCO, Oxygen of Carboxylic acid (protonated) SMARTS match = [CX3](=O)([OH1]) [OX2H1]([C](=O))
# [402, 404] = [[5], [1]]
bond 402 404 326.272386 1.36
```
* This type of comment is a match to the newer "amoeba21" database. 
* The SMARTS string match environment is given by [CX3](=O)([OH1]) [OX2H1]([C](=O)), where there is a space between the SMARTS for each atom.
* The first number in the bond parameter line is force constant and the second is the equilbrium bond length

##### Angle Parameter Definitions Example
```
# updated valence parameter database match, comments=O=C, Oxygen of carbonyl group, Acetic Acid C=O, sp2 carbon, carboxylic ester OCO, Oxygen of Carboxylic acid (protonated) SMARTS match = [OX1]=[CX3][OH1] [CX3](=O)([OH1]) [OX2H1]([C](=O))
# [406, 402, 404] = [[2], [5], [1]]
angle 406 402 404 109.848375 123.34
```
* The first number in the angle parameter line is force constant and the second is the equilbrium angle length

##### Stretch Bend Parameter Definitions Example
```
# updated valence parameter database match, comments=O=C, Oxygen of carbonyl group, Acetic Acid C=O, sp2 carbon, carboxylic ester OCO, Oxygen of Carboxylic acid (protonated) SMARTS match = [OX1]=[CX3][OH1] [CX3](=O)([OH1]) [OX2H1]([C](=O))
# [406, 402, 404] = [[2], [5], [1]]
strbnd 406 402 404 7.6289 7.6289
```


##### Out of Plane Bend Parameter Definitions Example
```
# updated valence parameter database match, comments=C=O, sp2 carbon, carboxylic ester OCO, Oxygen of Carboxylic acid (protonated) SMARTS match = [CX3](=O)([OH1]) [OX2H1]([C](=O))
# [404, 402] = [[1], [5]]
opbend 404 402 0 0 116.1422
```
* The first two class numbers are atoms in a trigonal center, the last two 0's are wild card atom classes for any other atom class in the trigonal center
* The last number is the opbend force constant

##### Torsion Parameter Definitions Example
```
# matching SMARTS from molecule  [['[*]~[*]~[*]~[*]', [1, 2, 3, 4]]] to SMARTS from parameter file [['[#6](-[H])(-[H])(-[H])-[#6](-[H])(-[H])(-[H])', [2, 1, 5, 6]]] with tinker type descriptions [[('H', '"Alkane H3C-"'), ('C', '"Alkane CH3-"'), ('C', '"Alkane CH3-"'), ('H', '"Alkane H3C-"')]]
# [403, 401, 402, 404] = [[3], [4], [5], [1]]
# Fitted from Fragment  SMARTS [#6](-[#6](-[#8]-[H])=[#8])(-[#7](-[H])-[H])(-[H])-[H] torsion atom indexes = 7,1,2,3 with smarts torsion indices 5,2,3,4 from fragment 5_1_Index_0.mol
# torsion % [#6](-[#6](-[#8]-[H])=[#8])(-[#7](-[H])-[H])(-[H])-[H] % 5,2,3,4 % -3.883,-0.434,4.077
torsion 403 401 402 404 -3.883 0.0 1 -0.434 180.0 2 4.077 0.0 3
```
* The line that starts with "Fitted from Fragment" , indicates which fragment the torsion parameters were derived from (from fragment 5_1_Index_0.mol) for debugging purposes
* torsion atom indexes = 7,1,2,3, indicates the atom indices that the torsion belongs too in the fragment molecule
* with smarts torsion indices 5,2,3,4 indicates the atom order in the SMARTS string corresponding to the torsion
* The torsion parameter line reads as "F Angle Number", where F is the force constant for the cosine term, Angle is the phase angle for the cosine term and Number is the number corresponding to which cosine term (can be up to 6).

##### Solute Parameter Definitions Example
```
#SOLUTE-SMARTS 408 [#1]([OH1](C=O))
SOLUTE 408 2.574 2.758 2.9054
```

##### Polarize Parameter Definitions Example
```
# updated valence parameter database match, comments=O on carbonyl group SMARTS match = [OX1]=[CX3]
# [406] = [[2]]
polarize           406          0.9138     0.3900 402
```

##### Multipole Parameter Definitions Example
```
# [404] = [[1]]
multipole   404  408  402              -0.46637
                                        0.01789    0.00000    0.22745
                                       -0.04708
                                        0.00000   -0.49060
                                       -0.09766    0.00000    0.53768

```
* The first line contain the monopole charge
* The second line contain the dipole
* The last lines contain the quadrupole matrix

#### Poltype Log File
* Contains useful information on the status of the program and which current step in the flow diagram it is in
* Will print errors at the bottom of the log file
* Example poltype.log shown below.
```
Mon Apr  4 11:52:41 2022 Running on host: node74.bme.utexas.edu
Mon Apr  4 11:52:41 2022 Atom Type Classification
Mon Apr  4 11:52:41 2022 QM Geometry Optimization
Mon Apr  4 11:52:41 2022 Calling: psi4 water_3D-opt_1.psi4 water_3D-opt_1.log path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:52:45 2022 Normal termination: logfile=/home/bdw2292/PoltypeJobs/SymmetryWater/Temp/water_3D-opt_1.log path=/home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:52:45 2022 Searching Database
Mon Apr  4 11:52:46 2022 Gas Phase Single Point for GDMA
Mon Apr  4 11:52:46 2022 Calling: psi4 water_3D-dma.psi4 water_3D-dma.log path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:52:49 2022 Normal termination: logfile=/home/bdw2292/PoltypeJobs/SymmetryWater/Temp/water_3D-dma.log path=/home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:52:49 2022 Gaussian Distributed Multipole Analysis (GDMA)
Mon Apr  4 11:52:49 2022 Calling: /opt/gdma/gdma-2.3.3/bin/gdma < water_3D.gdmain > water_3D.gdmaout path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:52:50 2022 Local Frame Symmetry Detection
Mon Apr  4 11:52:50 2022 Define Polarization Groups and Polarization Parameters
Mon Apr  4 11:52:50 2022 Calling: poledit 1 water_3D.gdmaout /home/bdw2292/poltype2/ParameterFiles/amoebabio18_header.prm < water_3D-peditin.txt path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:52:51 2022 Calling: potential 1 water_3D.xyz -k water_3D_prefitmultipole.key path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:52:52 2022 Gas Phase High Level Single Point for Multipole Refinement
Mon Apr  4 11:52:52 2022 Calling: psi4 water_3D-esp.psi4 water_3D-esp.log path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:20 2022 Normal termination: logfile=/home/bdw2292/PoltypeJobs/SymmetryWater/Temp/water_3D-esp.log path=/home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:20 2022 Calling: Generating CUBE File from PSI4
Mon Apr  4 11:53:20 2022 Calling: potential 2 water_3D_fortinker.cube path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:20 2022 Average Multipoles Via Symmetry
Mon Apr  4 11:53:20 2022 Calling: /home/bdw2292/poltype2/PoltypeModules/avgmpoles.pl water_3D_prefitmultipole.key water_3D.xyz water_3D-groups.txt water_3D.key_2 water_3D.xyz_2 401 path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:22 2022 Electrostatic Potential Optimization
Mon Apr  4 11:53:22 2022 Calling: potential 6 combined.xyz -k water_3D.key_2 combined.pot N 0.1 path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:23 2022 
Mon Apr  4 11:53:23 2022 =========================================================
Mon Apr  4 11:53:23 2022 Electrostatic Potential Comparison

Mon Apr  4 11:53:23 2022 Calling: potential 5 combined.xyz -k water_3D_postfitmultipole.key combined.pot N > RMSPD.txt path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:24 2022 RMSPD = 0.0877 Absolute tolerance is 1 kcal/mol and relative RMSPD=0.38% relative tolerance is 3%
Mon Apr  4 11:53:24 2022 Calling: minimize testbondangleequilvalues.xyz -k testbondangleequilvalues.key 0.1 > testbondangleequilvalues.out path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:24 2022 Calling: analyze testbondangleequilvalues.xyz_2 -k testbondangleequilvalues.key  d > testbondangleequilvaluesalz.out path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:24 2022 
Mon Apr  4 11:53:24 2022 =========================================================
Mon Apr  4 11:53:24 2022 Minimizing structure

Mon Apr  4 11:53:24 2022 Calling: minimize -k final.key final.xyz 0.1 > Minimized_final.out path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:25 2022 

Mon Apr  4 11:53:25 2022 =========================================================

Mon Apr  4 11:53:25 2022 Structure RMSD Comparison


Mon Apr  4 11:53:25 2022 Calling: superpose water_3D.xyz_2 final.xyz_2 1 N M N 0  > water_3D-superin.txt path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:25 2022 RMSD = 0.000000 Tolerance is 1
Mon Apr  4 11:53:25 2022 Calling: analyze water_3D.xyz_2 -k final.key em | grep -A11 Charge>MMDipole.txt path = /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:26 2022 Relative error of 0.0005361930294906766 for QMDipole 1.865 and 1.866 for MMDipole  tolerance = 0.5 /home/bdw2292/PoltypeJobs/SymmetryWater/Temp
Mon Apr  4 11:53:26 2022 Poltype Job Finished

```

#### OPENME Plots

<img src="README/Images/Torsion.PNG" width="50%">

*  Example plot of energy (left axis) and WBO (right axis) vs dihedral angle (bottom axis) torsion RMSE(red,blue)=.85 kcal/mol, relative RMSE=.42 kcal/mol. Yellow curve represents WBO, the green curve is AMOEBA prefit total energy, the blue curve is QM total energy, the pink curve is the fitting numerical spline + AMOEBA prefit total energy and the red curve is AMOEBA postfit total energy with an additional minimization with new parameters before evaluating the final AMOEBA energies.
*  The post fit red curve should be similar to the blue QM curve if fitting was successful (evaulated by RMSE). 
*  The pink curve and red curve should also be similar.

