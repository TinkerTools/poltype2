%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-062.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749597    0.104659    0.002363
 C   -0.666048   -0.014741   -0.000467
 H   -1.022542   -0.536980    0.900414
 H   -1.019709   -0.529535   -0.906735
 H   -1.120306    0.986946    0.002904
 H    1.128408   -0.788750   -0.000679

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

