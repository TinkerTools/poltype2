## Equilibriation Readme

## Equilibriation Protocol
* 1) For complexation, if restraint groups are not defined, automatically determine them.
* 2) For complexation, Compute center of mass between groups from minimized system, add keywords to keyfile.
* 3) For each temperature and restrain-position force constant (restrain position of protein), run MD in NVT. Default radius for restrain-position is 2 angstroms.
* 4) After reaching final temperature, run NPT dynamics. 
* 5) Average the box size from NPT dynamics output. 
* 6) Add new box sixe to keyfile (*_config.key)
* 7) If using group restraints during production MD, recompute COM average from .arc file and adjust keywords in file, else, compute the ideal distance, angle and torsion angle for restraints during production MD and compute the analytical correction from ensemble of distances, angles from equilbriation. 
* 8) Extract the last frame from XYZ file to be used for production dynamics.
