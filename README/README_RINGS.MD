## Rings README
<ol>
<li>	Find all non-aromatic ring torsions
<li>    For each set of non-aromatic torsions (separated by rings)
<ol>
    <li>Find all possible puckering locations </li>
<ol>
        <li>	Enumerate all torsions in ring of N-3 (how many torsions needed to span whole ring) </li>
        <li>	For each combination, if bond order >1, don’t use combination of torsions. If other torsions in combination are not neighbors of first torsion, then don’t use this combination. We want only consecutive torsions around ring to minimize restraints that could fight with each other. </li>
<ol>
    <li>Just use first combination. </li>
    <li>Update list of torsions used as input for fragmenter and torsiongenerator modules. </li>
    <li>Update list of torsions used as variable torsions (no torsion restraints added) during minimization for MM or QM. </li>
    <li>For each torsion in torsionset (group of torsions that will pucker for non-aromatic ring), assume a grid of 50 degrees with 10-degree increments. So effectively grid of 5x5x5… 5^(N-3) grid points for torsion set. </li>
    <li>Points from N-dimensional grid are uniformaly removed to keep a default of about ~40 QM points </li>
    <li>Update maxrange for each torsion (used to decide the angle increment every 12 degrees etc… ) in torsion set </li>
    <li>Compute the number of torsion parameters needed based on max range, number of torsions and number of cosine terms. </li>
    <li>Compute the number of datapoints needed based on number of parameters needed, then compute angle increment based on number of datapoints needed. </li>
    <li>For each torsion, use an initial guess (most torsions being fit do not use initial guess since they start as 0) of whatever is transferred from database of torsion parameters. </li>
<ol>
<ol>
