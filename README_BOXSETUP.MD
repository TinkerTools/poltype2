## PBC Box Setup

### Box Setup Protocol
* 1 Compute the cubic box size
* Take the longest dimension of protein or ligand (max (2*largest_cutoff+2*3, longest dimension+20A))
* 2 Compute the total number of waters based on box size
* 3 Compute the total number of physilogical ions based on listofsaltcons input
* 4 Compute the total number of neutralizing ions based on input system charge 
* 5 Compute the total atom number
* 6 Create solvent box based on total water number (using Tinker xyzedit)
* 7 Soak molecule in solvent box
* 8 Add physiological and neutralizing ions to solvent box
