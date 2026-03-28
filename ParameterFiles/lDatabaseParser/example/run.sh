# Assign charge flux parameters
python ../lAssignAMOEBAplusPRM.py -xyz phenol.xyz -key phenol.key -potent CF
# Assign polarizability parameters
python ../lAssignAMOEBAplusPRM.py -xyz phenol.xyz -key phenol.key -potent polar
# Assign bonded parameters, donot touch equilibrium b0, theta0
python ../lAssignAMOEBAplusPRM.py -xyz phenol.xyz -key phenol.key -potent bonded  
# Assign bonded parameters, with b0, theta0 being assigned as well
python ../lAssignAMOEBAplusPRM.py -xyz phenol.xyz -key phenol.key -potent bonded -konly no 
# Assign vdw parameters for amoeba ff
python ../lAssignAMOEBAplusPRM.py -xyz Me.xyz -key Me.prm -potent VDW