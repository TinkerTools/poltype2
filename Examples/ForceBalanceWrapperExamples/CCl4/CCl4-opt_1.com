%RWF=/scratch/bdw2292/Gau-CCl4/,100GB
%Nosave
%Chk=CCl4-opt_1.chk
%Mem=15GB
%Nproc=4
#P opt=(Cartesian,maxcycles=400) MP2/6-31G* MaxDisk=100GB 

CCl4 Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
Cl   -1.163000   -0.411800   -1.285400
Cl   -0.610600   -0.592500    1.565500
Cl    0.203200    1.768600    0.072700
Cl    1.570400   -0.764300   -0.352800
 C   -0.000000   -0.000000   -0.000000

