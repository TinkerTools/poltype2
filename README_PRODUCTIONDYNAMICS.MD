## Production Dynamics Readme

### Production Dynamics Protocol
* 1) For each electrostatic lambda, vdw lambda, restraint lambda, modify keywords in the keyfile and generate folders specific to each tuple of lambdas. 
* 2) If need to add ions to neutralize system, add multipole frame definitions for each ion in keyfile. Modify the charge according to the electrostatic lambda.
* 3) Copy the last frame from equilibriation to each folder.
* 4) Run production dynamics in NVT ensemble. 

