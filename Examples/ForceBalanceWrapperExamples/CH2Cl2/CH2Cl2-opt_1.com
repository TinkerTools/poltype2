%RWF=/scratch/bdw2292/Gau-CH2Cl2/,100GB
%Nosave
%Chk=CH2Cl2-opt_1.chk
%Mem=15GB
%Nproc=4
#P opt=(Cartesian,maxcycles=400) MP2/6-31G* MaxDisk=100GB 

CH2Cl2 Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
Cl   -1.457500    0.878700   -0.049100
Cl    1.486100    0.849400   -0.027200
 C   -0.000800   -0.139500   -0.009700
 H   -0.030700   -0.724500    0.935800
 H    0.003000   -0.863900   -0.849800

