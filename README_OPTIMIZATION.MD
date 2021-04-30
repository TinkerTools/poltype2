## Optimization Readme
<ol>
<li>	Create Psi4 OPT input file
<ol>
    <li>	Take input of coordinates, and total molecule charge, max opt cycles, QM method, QM basis set, scratch directory, number of threads, RAM. </li>
    <li>	Optional input of restrained bonds and angles (from vdw restrained QM optimization). </li>
    <li>	Use GAU LOOSE criteria for all optimizations. </li>
    <li>	If Iodine in molecule, then read input basis sets from file instead of keyword. </li>
    <li>	PCM is optional via keyword. </li>
</ol>
<li>	Create Gaussian OPT input .com file. </li>
<ol>
    <li>	Take input of coordinates, and total molecule charge, max opt cycles, QM method, QM basis set, scratch directory, number of threads, RAM. </li>
    <li>	Optional input of restrained bonds and angles (from vdw restrained QM optimization). </li>
    <li>	Use GAU LOOSE criteria for all optimizations. </li>
    <li>	If Iodine in molecule, then read input basis sets from file instead of keyword. </li>
    <li>	PCM is optional via keyword. </li>
    <li>	ModRedundant is default unless keyword specifies not to. </li>
<ol>
    <li>	In vdw optimization, don’t use mod redundant for dimers causes issues. </li>
</ol>
</ol>
<li>	Geometry optimization. </li>
<ol>
    <li>	If using Gaussian create Gaussian input file, if using Psi4 create Psi4 input file. </li>
    <li>	If external daemon is specified, then submit jobs via daemon else submit on local host. </li>
    <li>	Check for normal termination signal, if not found and error detected than crash program. </li>
    <li>	If using Gaussian call formchk program to generate fchk file. </li>
    <li>	If normal termination detected (number of cycles is >=maxoptcycles) then grab the final XYZ structure and output to -opt.xyz file. </li>
</ol>
<li>	Compare bond lengths. </li> 
<ol>
    <li>	If QM structure and input structure differ by tolerance of .2 angstroms for any bond then return signal to restart QM optimization with more opt cycles. </li>
</ol>
<li>	Check bond connectivity between final QM mol and input mol, if bond connectivity changed (such as hydrogen popping off etc…) then raise error and crash program. </li>
<li>	Final MM structure minimization after final parameters for comparison to QM structure. </li>
<ol>
    <li>	Check RMSD of final MM and QM structure in poltype and ensure value is less than 1. </li> 
</ol>
</ol>
