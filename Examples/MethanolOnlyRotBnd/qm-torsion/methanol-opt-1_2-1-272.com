%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-272.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749303    0.107916    0.005347
 C   -0.666711   -0.010908    0.000660
 H   -1.020163   -0.529582    0.904426
 H   -1.012847   -0.531857   -0.904416
 H   -1.132346    0.986859   -0.002920
 H    0.986164    0.928872    0.465503

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

