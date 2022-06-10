%RWF=/scratch/bdw2292/Gau-39_Index_0-farH_in17/,16.666666666666668GB
%Nosave
%Chk=39_Index_0-opt_1.chk
%Mem=8.333333333333334GB
%Nproc=1
#P opt=(Cartesian,maxcycles=400) MP2/6-31G* MaxDisk=16.666666666666668GB 

39_Index_0 Gaussian OPT Calculation on node44.bme.utexas.edu

0 1
 O    1.178200   -0.131000    0.062300
 C   -0.025500   -0.025300   -0.000600
 H   -0.422900    0.996200   -0.046600
 H   -0.729800   -0.840000   -0.015100

