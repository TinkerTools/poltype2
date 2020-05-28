## Torsion README

* Below is a high level flow of how torsion scanning works

* 2D torsion scanning is not yet supported.

* If doing a rigid rotation scan (all torsions around bond move together), all torsions around other bonds are restrained during optimization

* If doing a vibrational torsion scan, (slightly move one torsion around bond without moving others), all torsions around other bonds are restrained during optimization as well as all torsion besides the one being moved around rotatable bond

* Vibrational torsion scans will help partition the total QM-MM energy around a rotatable bond into UNIQUE parameter solutions for each individual torsion, giving a more accurate vibrational potential energy surface. This feature is not yet implmented in Poltype.

![Torsion Protocol](Images/TorsionProtocol.PNG)

* Below is more detailed flow of how files are generated

![Torsion Protocol](Images/TorsionProtocolDetailed.PNG)

* Below is a high level flow of how torsion fitting works

![Torsion Fitting Protocol](Images/TorsionFittingProtocol.PNG)

* Below is a plot of QM-MM and the coresponding fit vs dihedral angle

![QM-MM](Images/QM-MM.png)

* Below is a plot of total QM, MM energies vs dihedral angle

![Tor Energy](Images/TorEnergy.png)


