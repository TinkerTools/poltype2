%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-152.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749299    0.107426    0.001160
 C   -0.666494   -0.013365   -0.002473
 H   -1.016665   -0.538918    0.898207
 H   -1.025511   -0.526540   -0.908080
 H   -1.123753    0.987293    0.004286
 H    1.055524    0.139804   -0.919100

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

