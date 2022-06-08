%RWF=/scratch/bdw2292/Gau-CHCl3/,100GB
%Nosave
%Chk=CHCl3-opt_1.chk
%Mem=15GB
%Nproc=4
#P opt=(Cartesian,maxcycles=400) MP2/6-31G* MaxDisk=100GB 

CHCl3 Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
Cl   -0.817300    1.483200   -0.447200
Cl    1.658800    0.011800   -0.467800
Cl   -0.825900   -1.460200   -0.487700
 C   -0.022500   -0.001900    0.142400
 H    0.006900   -0.033000    1.260300

