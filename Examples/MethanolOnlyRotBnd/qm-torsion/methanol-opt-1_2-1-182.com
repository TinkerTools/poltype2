%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-182.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749213    0.107442    0.001648
 C   -0.666538   -0.012122   -0.001617
 H   -1.016984   -0.538736    0.898325
 H   -1.020230   -0.524531   -0.909270
 H   -1.127053    0.987435    0.006316
 H    1.014093    0.598611   -0.792302

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

