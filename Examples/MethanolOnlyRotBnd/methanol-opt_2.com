%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt_2.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=3,Loose) MP2/6-31G* MaxDisk=200GB 

methanol Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
 O    0.746152    0.122152   -0.000001
 C   -0.659980   -0.019687    0.000000
 H   -1.030589   -0.543576    0.893053
 H   -1.030588   -0.543617   -0.893032
 H   -1.078782    0.989156   -0.000018
 H    1.130627   -0.761056    0.000009

