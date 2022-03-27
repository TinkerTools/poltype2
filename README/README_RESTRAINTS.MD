## Restraints Readme

* Protein-ligand restraints
Protein-ligand restraints are usually used in bindng free energy calculations so that the ligand does not move away from the binding pocket even when it is decoupled from the system (vdw/ele-lambda = 0/0).
There is a missing step from the standard state (ligand that can freely move at the standard concentration) to the restrained state. The free energy change for this step can be evaluated either analytically or numerically.
This free energy change is often called DeltaG_correction or DeltaG_restraint.
DeltaG_restraint should be positive for binding/solvation, and negative for unbinding/desolvation.
* In output CSV file the value is referenced under ΔGᵃⁿᵃᶜᵒᵐᵖᶜᵒʳʳ 


* Types of restraint
group-restraint: one distance restraint between a group of atoms in the ligand and a group of atoms in the protein. This is used during equilbriation. 
Boresch restraint: 6 restraints defined by 3 atoms in the ligand and 3 atoms in the protein, which restrains the relative translation and rotation. This is used by default during production dynamics.
# https://doi.org/10.1021/jp0217839


For both restraints, DeltaG_restraint for the recommended/default restraint keyword is calculated.
For group retraint, the user needs to uncomment the "restrain-groups" keyword or specify custom force constant and r0.

