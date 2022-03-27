## Vdw README
<ol>
<li>	Convert cartesian probe dimer XYZ to Tinker XYZ </li>
<li>	Generate initial probe structures </li>
<ol>
    <li>	Start with dictionary of element to vdw radii </li>
    <li>	Grab symmetry types of input molecule </li>
    <li>	Grab symmetry types of probe molecule </li>
<ol>
        <li>	Right now, either water or just molecule again (homo dimer) </li>
</ol> 
    <li>	Grab symmetry types of probe molecule </li>
    <li>	Grab only one instance of atom indices corresponding to molecule/probe atom types </li>
    <li>	For every molecule atom-probe atom </li>
<ol>
        <li>	Excluding SP3 -SP3 (like in homodimer), excluding H-H atoms </li>
        <li>	Generate initial dimer </li>
</ol>
</ol>
<li>	Generate initial dimer </li>
<ol>
    <li>	Generate reference distances (angstroms) </li>
<ol>
        <li>	For target probe atom- target molecule atom, use combining rule for vdw from each elements vdw radius and use that as target distance </li>
        <li>	For every other probe atom- molecule atom, treat target distance as >=1 </li>
        <li>	Upper and lower bounds for each distance, for steric, trivially high bound, for target upper and lower bound is target distance </li>
</ol>
    <li>	Generate reference angles </li>
<ol>
        <li>	Three types of triangles between probe and molecule can be made to put extra constraints on cost function (orientation of probe relative to dimer). Acceptor atom is molecule target atom index, donor atom is probe target atom index </li>
<ol>
<li>	A=Neighbor of donor atom, B= donor atom, C= acceptor atom </li>
<ol>
    <li>	Special case for water, let O point towards molecule and the two H point away from molecule </li>
</ol>
<li>	A=donor atom, neighbor of acceptor atom, acceptor atom </li>
<ol>
    <li>	If multiple acceptor neighbors and hybridization of acceptor is SP2, then set angle of this triangle to 90 degrees (think for hydrogen pointing towards SP2 atom perpendicular instead of at angle) </li>
</ol>
<li>	A=Neighbor of Neighbor of acceptor atom, B= Neighbor of acceptor atom, C= donor atom </li>
<ol>
    <li>	If only one acceptor neighbor and multiple neighbors of acceptor (think SP2 N with hydrogen as acceptor, and oxygen in water as donor, then this will keep the water probe in same plane as SP2 N) </li>
</ol>
</ol>
</ol>
    <li>	Search on 3D mesh grid for optimal probe location near target atom on molecule while avoiding sterics with molecule. </li>
    <li>	Convert final cartesian XYZ probe to tinker XYZ probe </li>
    <li>	Minimize tinker XYZ with the above distance and angle restraints with some tolerance </li>
<ol>
        <li>	    distanceforceconstant=5 # kcal/mol/ang^2 </li>
        <li>	    angleforceconstant=.1 # kcal/mol/deg^2 </li>
        <li>	    Generate lower and upper bounds for input distances and angles </li>
<ol>
<li>	+/- 10% of distance restraint </li>
<li>	+/- 10 degrees </li>
</ol>
</ol>
</ol>
<li>	For each tinker minimized initial probe structure </li>
<ol>
    <li>	Try geometry optimization (by default dont do QM OPT, just use MM OPT) </li>
<ol>
        <li>	If this fails, then remove this probe atom – molecule atom pair from fitting </li>
        <li>	If successful grab output XYZ and then convert to tinker XYZ </li>
</ol>
    <li>	Then move dimer about the QM minima (to sample many points along QM vs distance surface) </li>
    <li>	Then generate input SP files for Gaussian/Psi4 </li>
    <li>	Execute Gaussian/Psi4 SP jobs </li>
</ol>
<li>	Move dimer about MM (or QM if want QM minima) minima </li>
<ol>
    <li>	Default array of .8,.9,1,1.1,1.2 </li>
    <li>	The initial distance between target atom on the molecule and the probe atom is scaled to each value in array by shifting the probe </li>
    <li>	Generate output tinker xyz file for dimer </li>
</ol>
<li>	Check if all .png images have been generated (if already tried doing the fits), if not then proceed to vdw fitting </li>
<ol>
    <li>	Combine molecule atom – probe atom pairs that need to be fit together, such as same type probed by more than one type </li>
    <li>	Read counterpose QM data and write to QM.dat file </li>
    <li>	Grab initial vdw parameters from key_5 file and determine bounds </li>
<ol>
        <li>	Radius is +/- 10% (radius) </li>
        <li>	Depth is +/- 10% (depth) </li>
        <li>	Reduction is -.2 (absolute value) and upper bound of 1 </li>
</ol>
    <li>	Write initial guess and bounds to INITIAL.PRM file </li>
    <li>	Optimize VDW parameters </li>
</ol>
<li>	Optimize VDW Parameters </li>
<ol>
    <li>	Read from INITIAL.PRM file </li>
<ol>
        <li>	Read initial guess radius, depth and reduction factor </li>
        <li>	Reduction factor is 1 if not an atom type with a valence of 1 </li>
        <li>	Read min/max bounds for radius, depth and reduction factor </li>
        <li>	Put into format for optimizer </li>
<ol>
        <li>	Using least squares optimizer from scipy </li>
</ol>
</ol>
</ol>

<li>	Vdw cost function </li>
<ol>
    <li>	Input vdw parameters </li>
    <li>	Read QM data energies </li>
<ol>
        <li>	Normalize </li>
        <li>	Screen energies > 15 kcal/mol </li>
</ol>
    <li>	Write parameters to key file </li>
    <li>	Run analyze </li>
    <li>	Read MM data energies
<ol>
        <li>	Normalize </li>
        <li>	Screen energies > 15 kcal/mol </li>
</ol>
    <li>	Return weightlist*(QM-MM) </li>
<ol>
        <li>	If first time fitting, then weights are all 1 </li>
        <li>	Else use Boltzmann equation on QM data with default RT of 8 </li>
</ol>
</ol>
<li>	Plot Energy Vs Distance for each probe atom-molecule atom pair </li>
<ol>
    <li>	Plot final Parameter values </li>
    <li>	Plot RMSE </li>
</ol>
<li>	Plot QM Vs MM Energy </li>
<ol>
    <li>	Plot best fit line </li>
    <li>	Plot r^2 value </li>
    <li>	Plot RMSE value </li>
    <li>	If RMSE is too high after Boltzmann fitting, then crash program default of 1.6 RMSE and relative RMSE of .2 </li>
</ol>
</ol>
