%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-302.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749493    0.107064    0.005035
 C   -0.666231   -0.012340    0.001485
 H   -1.023828   -0.531796    0.903581
 H   -1.012886   -0.531828   -0.904064
 H   -1.126613    0.987331   -0.000539
 H    1.011066    0.591969    0.803903

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

