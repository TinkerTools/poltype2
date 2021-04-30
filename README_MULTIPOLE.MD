## Multipoles

![Multipole Frames](Images/MpoleFrames.PNG)


* Top left: bisector
* Top middle: z-then-x
* Top right z-only
* Bottom left: bisector-then-x
* Bottom middle: z-then-x-then-y
* Bottom right: trisector

* examples of these coordinate systems can be found in README_MANIFEST, see SymmetryExamples

### Multipole Parameters
```
M1 M2 M3      q
              Dx  Dy  Dz
              Qxx
              Qxy  Qyy
              Qxz  Qyz  Qzz
```

### Z-then-X
* M1 M2 M3 

### bisector
* M1 M2 -M2, or M1 -M2 M2, or M1 -M2 -M2

### bisector-then-X
* M1 M2 -M3 -M3

* typically Dy, Qxy, Qyz components are zero, unless it is chiral such as C_alpha in protein
## Multipole.py Protocol

<ol>
<li>	Generate poledit input file </li>
<ol>
   <li> Let poledit handle frame definitions </li>
   <li>	Let poledit handle polarizabilities except if aromatic carbon and aromatic hydrogen </li>

   <li>	Tell poledit to make polariable groups via rotatable bonds </li>
<ol>
	<li> Start with all atoms same group </li>
	<li> Separate into rotatable bonds </li>
	<li> Remove bonds involving univalent atoms </li>
	<li> Remove bonds internal to aromatic ring </li>
	<li> Remove bonds with SP carbon </li>
	<li> Remove C=C bond of terminal alkene </li>
	<li> Remove C-O bonds of alcohol, ether </li>
	<li> Remove C-C in ethane </li>
	<li> Remove C-O bonds of carboxylic acid , ester </li>
	<li> Remove C-N bonds of alkyl amine </li>
	<li> Remove C-N bonds of amide, urea, amidine, guanidinium </li>
	<li> Remove any P-X, S-X or with X = (N or O) </li>
</ol>

   <li>	Tell poledit not to condense equivalent types (we will do this ourselves with symmetry classes) </li>

   <li>	Tell poledit to zero out components in dipole and quadrupole due to symmetry in frame </li>
</ol>

<li>    Generate input for GDMA </li>
<ol>
   <li>	Set radius of integrating electron potential for elements H, S, P, Cl, Br, and I. Tell GDMA to use protocol from first GDMA version. </li>

</ol>

<li>    Execute GDMA </li>
<ol>
   <li>	If using Gaussian then use executable else take output from Psi4 DMA SP log file </li>

<li>    Prepend key file </li>

   <li>	Turn off parameters for everything except polarize and multipole (for esp fitting) </li>

</ol>

<li>    Average multipoles </li>

   <li>	Generate inputs for perl script (.xyz, .key and all atom indices belonging to same symmetry type) </li>

   <li>	Call external perl script to parse and average all values in multipoles via same symmetry types. Generate .xyz_2 and .key_2 files. </li>

</ol>
</ol>


## A Special example of Nitrate

Note NO3- is fairly flat and rigid so we are able to use Z-then-X for both N and O. For N, all dipole components have to be 0 and the Qxx and Qzz are equal (x and Z are in the plane, y is normal to the plane). For O, the "2 1 2" also works because the direction of "X" or "Y" is not affecting the (sign) of Qxx or Qyy.
```
--------AMOEBA Parameters for Nitrate Ion-----


atom 1 1 N "Nitrate N" 7 14.007 3

atom 2 2 O "Nitrate O" 8 15.999 1

vdw 1 3.7100 0.1100

vdw 2 3.5100 0.1120

bond 1 2 390.0 1.2606

angle 2 1 2 155.00 120.00

strbnd 2 1 2 18.70 18.70

opbend 2 1 2 2 172.00

multipole 1 2 2 1.07363

                  0.00000 0.00000 0.00000
                  
                  0.14407
                  
                  0.00000 -0.28814
                  
                  0.00000 0.00000 0.14407

multipole 2 1 2 -0.69121

                 0.00000 0.00000 0.04941
                 
                -0.50067
                
                 0.00000 -0.08271
                 
                 0.00000 0.00000 0.58338

polarize 1 1.0730 0.3900 2

polarize 2 0.8370 0.3900 1
```
