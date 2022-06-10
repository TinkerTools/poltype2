%RWF=/scratch/bdw2292/Gau-74_Index_0-farH_in17/,16.666666666666668GB
%Nosave
%Chk=74_Index_0-opt_1.chk
%Mem=8.333333333333334GB
%Nproc=1
#P opt=(Cartesian,maxcycles=400) MP2/6-31G* MaxDisk=16.666666666666668GB 

74_Index_0 Gaussian OPT Calculation on node44.bme.utexas.edu

1 1
 H   -0.650500   -0.680300   -0.451600
 N    0.015600   -0.005300    0.003400
 H   -0.516800    0.445700    0.785200
 H    0.283300    0.754100   -0.664000
 H    0.868300   -0.514200    0.327000

