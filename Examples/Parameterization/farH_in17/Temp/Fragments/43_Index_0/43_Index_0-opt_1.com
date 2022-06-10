%RWF=/scratch/bdw2292/Gau-43_Index_0-farH_in17/,16.666666666666668GB
%Nosave
%Chk=43_Index_0-opt_1.chk
%Mem=8.333333333333334GB
%Nproc=1
#P opt=(Cartesian,maxcycles=400) MP2/6-31G* MaxDisk=16.666666666666668GB 

43_Index_0 Gaussian OPT Calculation on node44.bme.utexas.edu

0 1
 H   -0.857800   -0.502100   -0.085000
 N    0.006700   -0.006000    0.259300
 H    0.888800   -0.464000   -0.089300
 H   -0.037700    0.972100   -0.084900

