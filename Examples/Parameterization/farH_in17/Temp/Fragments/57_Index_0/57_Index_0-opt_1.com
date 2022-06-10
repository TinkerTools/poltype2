%RWF=/scratch/bdw2292/Gau-57_Index_0-farH_in17/,16.666666666666668GB
%Nosave
%Chk=57_Index_0-opt_1.chk
%Mem=8.333333333333334GB
%Nproc=1
#P opt=(Cartesian,maxcycles=400) MP2/6-31G* MaxDisk=16.666666666666668GB 

57_Index_0 Gaussian OPT Calculation on node44.bme.utexas.edu

0 1
 H   -0.742200   -0.246700    0.758900
 C    0.026400    0.011200   -0.012600
 H   -0.362500    0.908400   -0.531200
 H    0.069600   -0.859900   -0.700600
 H    1.008700    0.186900    0.485600

