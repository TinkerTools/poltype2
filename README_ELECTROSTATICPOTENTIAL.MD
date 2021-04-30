## Electrostatic Potential Readme
<ol>
<li>	Create Psi4 Input file used for DMA SP computation </li>
<ol>
    <li>Take inputs of charge, molecule coordinates, RAM, scratch directory, number of threads, QM method, QM basisset </li>
    <li>If Iodine is in molecule, then instead of default QM basiset, read from input files containing Iodine ECP parameters </li>
<ol>
        <li>	Note Psi4 has issues with ECP last time checked. So poltype defaults to gaussian for everything if gaussian is detected in path and also Iodine is in molecular formula.</li>
</ol>
    <li>Frozen electron core is on by default </li>
    <li>COM frame of reference is on by default (compute properties from COM) </li>
    <li>Generate esp cube file for visualizing in GaussView / VMD </li>
    <li>Compute the dipole vector, ESP grid, WBO index </li>
    <li>Compute fchk file </li>
    <li>GDMA via psi4 after DMA SP finished </li>
    <li>Delete scratch files after words </li>
</ol>
<li>	Create Psi4 Input file used for ESP SP computation </li>
<ol>
    <li>Same as function for DMA input file, but different inputs for QM basis set/QM method </li>
    <li>Also don’t compute GDMA via psi4 after DMA SP finished </li>
</ol>
<li>	Generate COM file for Gaussian DMA/ESP SP computations </li>
<ol>
    <li>Take inputs of charge, molecule coordinates, RAM, scratch directory, number of threads, QM method, QM basis set. </li>
    <li>If Iodine is in molecule, then instead of default QM basis set, read from input files containing Iodine ECP parameters. </li>
</ol>
<li>	Call and compute DMA SP, ESP SP </li> 
<ol>
    <li>If using gaussian, then generate input COM file, otherwise generate input Psi4 file. </li>
    <li>If user specified external daemon, then use that to submit jobs, otherwise call QM job on local host. </li>
    <li>Check for termination signal after calling psi4/gaussian, if program did not terminate and has error than raise error in poltype and crash program. </li>
</ol>
<li>	Generate ESP .pot file for tinker prior to ESP multipole fitting </li>
<ol>
    <li>If used psi4 to generate DMA and ESP QM data, then read from psi4’s grid_esp.dat file and add contents to a cube file (only what tinker needs not compatible as an actually gaussian cube since missing information in file). Also read from tinkers .grid file and add to cube file. </li>
<ol>
    <li>	Otherwise use cubegen program that comes with gaussian on fchk file, using .grid file from first tinker potential program call (generates .key, .xyz and .grid) </li>
</ol>
    <li>Run tinkers’ potential program on cube file to extract information and create the .pot file. </li>
</ol>
<li>	Electrostatic Potential Fitting </li>
<ol>
    <li>Call tinkers potential program with .pot , .xyz_2 and .key_2 as inputs. </li>
</ol>
<li>	Compare final multipole parameters to QM potential grid file </li>
<ol>
    <li>Compare via tinkers potential program. </li>
    <li>Check Root Mean Square Potential Difference between final multipole fitted parameters and QM potential grid. Raise error if greater than tolerance of 1. Also must fail relative RMSPD criteria of RMSPD/atom number with tolerance of .013 </li>
</ol>
<li>	Grab QM dipoles </li> 
<ol>
    <li>From Psi4 ESP log file or Gaussian ESP log file. </li>
    <li>Gaussian dipole vector is converted to COM frame dipole vector (tinker computes in COM frame). </li>
<ol>
	<li>Compute COM position vector and COC (center of charge) position vector, then compute displacement vector needed to convert the dipole vector </li>
</ol>
</ol>
<li>	Grab MM dipoles </li>
<ol>
    <li>Call analyze tinker program. </li>
</ol>
<li>	Compare QM/MM dipoles and if they differ more than a tolerance of .5, then crash poltype and raise error </li>
</ol>
