%RWF=/scratch/bdw2292/Gau-water_3D/,317GB
%Nosave
%Chk=water_3D-opt_1.chk
%Mem=49GB
%Nproc=12
#P opt=(Cartesian,maxcycles=400) MP2/6-31G* MaxDisk=317GB 

water_3D Gaussian OPT Calculation on node74.bme.utexas.edu

0 1
 O    0.009400    0.397700    0.000000
 H    0.758600   -0.216900    0.000000
 H   -0.768000   -0.180800    0.000000

