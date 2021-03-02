%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-332.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749546    0.106190    0.004767
 C   -0.666255   -0.013798    0.002041
 H   -1.028870   -0.533895    0.902198
 H   -1.013624   -0.531946   -0.904034
 H   -1.122520    0.987369    0.001447
 H    1.051823    0.131080    0.926581

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

